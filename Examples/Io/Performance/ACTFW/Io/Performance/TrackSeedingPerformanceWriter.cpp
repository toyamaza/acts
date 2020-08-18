// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Performance/TrackSeedingPerformanceWriter.hpp"

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SpacePoint.hpp"
#include <Acts/Utilities/Helpers.hpp>

#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using HitParticlesMap = FW::IndexMultimap<ActsFatras::Barcode>;

FW::TrackSeedingPerformanceWriter::TrackSeedingPerformanceWriter(
    FW::TrackSeedingPerformanceWriter::Config cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputSeeds, "TrackSeedingPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_effPlotTool(m_cfg.effPlotToolConfig, lvl) {
  // Input track and truth collection name
  if (m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument("Missing input seeds collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
  if (m_cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path = joinPaths(m_cfg.outputDir, m_cfg.outputFilename);
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (not m_outputFile) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }

  // initialize the residual and efficiency plots tool
  m_effPlotTool.book(m_effPlotCache);
}

FW::TrackSeedingPerformanceWriter::~TrackSeedingPerformanceWriter() {
  m_effPlotTool.clear(m_effPlotCache);

  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode FW::TrackSeedingPerformanceWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_effPlotTool.write(m_effPlotCache);

    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

void FW::TrackSeedingPerformanceWriter::printSeed(
    const Acts::Seed<SpacePoint>* seed) const {
  const SpacePoint* sp = seed->sp()[0];
  std::cout << "Sp Id # " << sp->Id() << ": # of particles "
            << sp->particles.size() << ": ";
  std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
            << sp->z()
            << ") particle = " << sp->particles[0].particleId.particle()
            << ", barcode: " << sp->particles[0].particleId;
  sp = seed->sp()[1];
  std::cout << "; Sp Id # " << sp->Id() << ": # of particles "
            << sp->particles.size() << ": ";
  std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
            << sp->z()
            << ") particle = " << sp->particles[0].particleId.particle()
            << ", barcode: " << sp->particles[0].particleId;
  sp = seed->sp()[2];
  std::cout << "; Sp Id # " << sp->Id() << ": # of particles "
            << sp->particles.size() << ": ";
  std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
            << sp->z()
            << ") particle = " << sp->particles[0].particleId.particle()
            << ", barcode: " << sp->particles[0].particleId << " with count "
            << sp->particles[0].hitCount;
  std::cout << std::endl << std::endl;
}

size_t FW::TrackSeedingPerformanceWriter::seedNumParticles(
    const Acts::Seed<SpacePoint>* seed,
    std::set<ActsFatras::Barcode>& particlesFoundBySeeds,
    size_t& nDuplicateSeeds) const {
  const SpacePoint* sp0 = seed->sp()[0];
  const SpacePoint* sp1 = seed->sp()[1];
  const SpacePoint* sp2 = seed->sp()[2];
  std::set<ActsFatras::Barcode> particles0;
  std::set<ActsFatras::Barcode> particles1;
  std::set<ActsFatras::Barcode> particles2;
  for (size_t i = 0; i < sp0->particles.size(); i++) {
    particles0.insert(sp0->particles[i].particleId);  // insert particle barcode
  }
  for (size_t i = 0; i < sp1->particles.size(); i++) {
    particles1.insert(sp1->particles[i].particleId);
  }
  for (size_t i = 0; i < sp2->particles.size(); i++) {
    particles2.insert(sp2->particles[i].particleId);
  }
  std::set<ActsFatras::Barcode> tmp;
  set_intersection(particles0.begin(), particles0.end(), particles1.begin(),
                   particles1.end(), std::inserter(tmp, tmp.end()));
  std::set<ActsFatras::Barcode> particlesInCommon;
  set_intersection(particles2.begin(), particles2.end(), tmp.begin(), tmp.end(),
                   std::inserter(particlesInCommon, particlesInCommon.end()));

  for (const auto& particle : particlesInCommon) {
    if (particlesFoundBySeeds.find(particle) != particlesFoundBySeeds.end()) {
      // this particle was already found
      nDuplicateSeeds++;
    }
    particlesFoundBySeeds.insert(particle);
  }
  return particlesInCommon.size();
}

FW::ProcessCode FW::TrackSeedingPerformanceWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<std::vector<Acts::Seed<SpacePoint>>>& seedVector) {
  size_t nSeeds = 0;
  for (auto& regionVec : seedVector) {
    nSeeds += regionVec.size();
  }

  // Read truth particles from input collection
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);

  size_t nParticles = particles.size();
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Set that helps keep track of the number of duplicate seeds. i.e. if there
  // are three seeds for one particle, this is counted as two duplicate seeds.
  std::set<ActsFatras::Barcode> particlesFoundBySeeds;
  size_t nDuplicateSeeds = 0;  // number of duplicate seeds as defined above
  size_t nTrueSeeds = 0;       // true seed means it contains a particle
  size_t maxNumToPrint = 0;    // 0 means print nothing
  size_t printCounter = 0;     // for debuging purposes
  for (auto& regionVec : seedVector) {
    for (size_t i = 0; i < regionVec.size(); i++) {
      const Acts::Seed<SpacePoint>* seed = &regionVec[i];
      size_t prtsInCommon =
          seedNumParticles(seed, particlesFoundBySeeds, nDuplicateSeeds);
      if (prtsInCommon > 0) {
        nTrueSeeds++;
      }
    }
    // prints first few seeds for debugging purposes
    if (printCounter < maxNumToPrint) {
      for (size_t i = 0; i < regionVec.size(); i++) {
        if (printCounter >= maxNumToPrint) {
          break;
        }
        const Acts::Seed<SpacePoint>* seed = &regionVec[i];
        printSeed(seed);
        printCounter++;
      }
    }
  }
  ACTS_INFO("Number of seeds generated: " << nSeeds)
  ACTS_INFO("Number of true seeds generated: " << nTrueSeeds)
  ACTS_INFO("Number of duplicate seeds generated: " << nDuplicateSeeds)
  ACTS_INFO("Fake rate (nSeeds - nTrueSeeds) / nSeeds --- "
            << 100 * (nSeeds - nTrueSeeds) / nSeeds << "%")
  ACTS_INFO("Technical Effeciency (nTrueSeeds - nDuplicateSeeds / nSeeds) --- "
            << 100 * (nTrueSeeds - nDuplicateSeeds) / nSeeds << "%")
  ACTS_INFO("Duplicate rate (nDuplicateSeeds / nTrueSeeds) --- "
            << (100 * nDuplicateSeeds) / nSeeds << "%")
  ACTS_INFO("Effeciency: (nTrueSeeds - nDuplicateSeeds)/nParticles = "
            << nTrueSeeds - nDuplicateSeeds << " / " << nParticles << " = "
            << 100 * (nTrueSeeds - nDuplicateSeeds) / nParticles << "%")

  // Fill the effeciency plots
  for (const auto& particle : particles) {
    const auto it = particlesFoundBySeeds.find(particle.particleId());
    if (it != particlesFoundBySeeds.end()) {
      // when the trajectory is reconstructed
      // 0 particle id means the particle is a fake one.
      m_effPlotTool.fill(m_effPlotCache, particle, it->particle() != 0);
    } else {
      // when the trajectory is NOT reconstructed
      m_effPlotTool.fill(m_effPlotCache, particle, false);
    }
  }
  return ProcessCode::SUCCESS;
}