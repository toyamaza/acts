// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TestSeedAlgorithm.hpp"

#include "ACTFW/EventData/IndexContainers.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimIdentifier.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ACTFW/Validation/ProtoTrackClassification.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Seeding/ATLASCuts.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePoint.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <boost/type_erasure/any_cast.hpp>

using HitParticlesMap = FW::IndexMultimap<ActsFatras::Barcode>;

FW::TestSeedAlgorithm::TestSeedAlgorithm(
    const FW::TestSeedAlgorithm::Config& cfg, Acts::Logging::Level level)
    : FW::BareAlgorithm("TestSeedAlgorithm", level), m_cfg(std::move(cfg)) {
  // Is inputClusters already checked by base constructor?
  // I think this is only true for the writer
  if (m_cfg.inputClusters.empty()) {
    throw std::invalid_argument(
        "Missing clusters input collection with the hits");
  }
  if (cfg.inputHitParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
}

// Technically, space points can have multiple particles that are a part of
// them, so seedNumParticles loops over all possible combinations of particles
// in the seeds to find how many particles are in common.
// TODO: Replace code with set intersection so that it's much more simple.
// Returns the number particles that are a part of all
// 3 spacePoints in the seed. 0 means it's a fake seed.
std::size_t FW::TestSeedAlgorithm::seedNumParticles(
    const Acts::Seed<SpacePoint>* seed) const {
  std::size_t numTrueSeeds = 0;
  const SpacePoint* sp0 = seed->sp()[0];
  const SpacePoint* sp1 = seed->sp()[1];
  const SpacePoint* sp2 = seed->sp()[2];
  for (std::size_t a = 0; a < sp0->particles.size(); a++) {
    for (std::size_t b = 0; b < sp0->particles.size(); b++) {
      if (sp0->particles[a].particleId.particle() ==
          sp1->particles[b].particleId.particle()) {
        for (std::size_t c = 0; c < sp0->particles.size(); c++) {
          if (sp1->particles[b].particleId.particle() ==
              sp2->particles[c].particleId.particle()) {
            numTrueSeeds++;
          }
        }
      }
    }
  }
  return numTrueSeeds;
}

// returns a space point with a particle barcode stored in particles for each
// particle that made this space point
SpacePoint* FW::TestSeedAlgorithm::readSP(
    std::size_t hit_id, const Acts::GeometryID geoId,
    const Acts::PlanarModuleCluster& cluster,
    const HitParticlesMap& hitParticlesMap, const AlgorithmContext& ctx) const {
  const auto& parameters = cluster.parameters();
  Acts::Vector2D localPos(parameters[0], parameters[1]);
  Acts::Vector3D globalFakeMom(1, 1, 1);
  Acts::Vector3D globalPos(0, 0, 0);
  // transform local into global position information
  cluster.referenceObject().localToGlobal(ctx.geoContext, localPos,
                                          globalFakeMom, globalPos);
  float x, y, z, r, varianceR, varianceZ;
  x = globalPos.x();
  y = globalPos.y();
  z = globalPos.z();
  r = std::sqrt(x * x + y * y);
  varianceR = 0;  // initialized to 0 becuse they don't affect seeds generated
  varianceZ = 0;  // initialized to 0 becuse they don't affect seeds generated

  // get truth particles that are a part of this space point
  std::vector<FW::ParticleHitCount> particleHitCount;
  for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hit_id))) {
    auto particleId = hitParticle.second;
    // search for existing particle in the existing hit counts
    auto isSameParticle = [=](const ParticleHitCount& phc) {
      return (phc.particleId == particleId);
    };
    auto it = std::find_if(particleHitCount.begin(), particleHitCount.end(),
                           isSameParticle);
    // either increase count if we saw the particle before or add it
    if (it != particleHitCount.end()) {
      it->hitCount += 1;
    } else {
      particleHitCount.push_back({particleId, 1u});
    }
  }

  SpacePoint* sp = new SpacePoint{
      x, y, z, r, geoId.layer(), varianceR, varianceZ, particleHitCount};
  return sp;
}

FW::ProcessCode FW::TestSeedAlgorithm::execute(
    const AlgorithmContext& ctx /*,
    const FW::GeometryIdMultimap<Acts::PlanarModuleCluster>& clusters*/) const {
  // read in the hits
  const auto& clusters =
      ctx.eventStore.get<FW::GeometryIdMultimap<Acts::PlanarModuleCluster>>(
          m_cfg.inputClusters);
  // read in the map of hitId to particleId truth information
  const HitParticlesMap hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputHitParticlesMap);

  // create the space points
  std::size_t clustCounter = 0;
  std::size_t numIgnored = 0;
  std::vector<const SpacePoint*> spVec;
  // since clusters are ordered, we simply count the hit_id as we read
  // clusters. Hit_id isn't stored in a cluster. This is how
  // CsvPlanarClusterWriter did it.
  std::size_t hit_id = 0;
  for (const auto& entry : clusters) {
    Acts::GeometryID geoId = entry.first;
    const Acts::PlanarModuleCluster& cluster = entry.second;
    std::size_t volnum = geoId.volume();
    std::size_t laynum = geoId.layer();

    // only uses innermost 3 pixel barrel-like layers on volume 8.
    if (volnum == 8) {
      if (laynum == 2 || laynum == 4 || laynum == 6) {
        SpacePoint* sp = readSP(hit_id, geoId, cluster, hitParticlesMap, ctx);
        spVec.push_back(sp);
        clustCounter++;
      } else {
        numIgnored++;
      }
    } else {
      numIgnored++;
    }
    hit_id++;
  }

  Acts::SeedfinderConfig<SpacePoint> config;
  // silicon detector max
  config.rMax = 160.;
  config.deltaRMin = 5.;
  config.deltaRMax = 160.;
  config.collisionRegionMin = -250.;
  config.collisionRegionMax = 250.;
  config.zMin = -2800.;
  config.zMax = 2800.;
  config.maxSeedsPerSpM = 5;
  // 2.7 eta
  config.cotThetaMax = 7.40627;
  config.sigmaScattering = 1.00000;

  config.minPt = 500.;
  config.bFieldInZ = 0.00199724;

  config.beamPos = {-.5, -.5};
  config.impactMax = 10.;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  Acts::SeedFilterConfig sfconf;
  Acts::ATLASCuts<SpacePoint> atlasCuts = Acts::ATLASCuts<SpacePoint>();
  config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
      Acts::SeedFilter<SpacePoint>(sfconf, &atlasCuts));
  Acts::Seedfinder<SpacePoint> a(config);

  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SpacePoint& sp, float, float, float) -> Acts::Vector2D {
    return {sp.varianceR, sp.varianceZ};
  };

  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;

  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid =
      Acts::SpacePointGridCreator::createGrid<SpacePoint>(gridConf);
  auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(), spVec.end(), ct,
                                                 bottomBinFinder, topBinFinder,
                                                 std::move(grid), config);

  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector;
  auto start = std::chrono::system_clock::now();
  auto groupIt = spGroup.begin();
  auto endOfGroups = spGroup.end();

  // actually executues the seed finding algoirthm here
  for (; !(groupIt == endOfGroups); ++groupIt) {
    seedVector.push_back(a.createSeedsForGroup(
        groupIt.bottom(), groupIt.middle(), groupIt.top()));
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "time to create seeds: " << elapsed_seconds.count() << std::endl;
  std::cout << "Number of regions: " << seedVector.size() << std::endl;
  std::cout << "Now we start counting seeds" << std::endl;
  int numSeeds = 0;
  std::size_t emptySeedCount = 0;
  for (auto& regionVec : seedVector) {
    numSeeds += regionVec.size();
    if (regionVec.size() == 0) {
      emptySeedCount++;
    }
  }
  std::size_t maxNumToPrint = 10;  // for debuging purposes
  std::size_t printCounter = 0;    // for debuging purposes
  std::size_t printCounter2 = 0;   // for debuging purposes
  std::size_t numTrueSeeds = 0;
  for (auto& regionVec : seedVector) {
    for (size_t i = 0; i < regionVec.size(); i++) {
      const Acts::Seed<SpacePoint>* seed = &regionVec[i];
      if (seedNumParticles(seed) > 0) {
        numTrueSeeds++;
        if (printCounter2 < maxNumToPrint) {
          std::cout << "We have a true seed: ";
          const SpacePoint* sp = seed->sp()[0];
          std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                    << sp->z() << ") hit_id? = " << sp->particles[0].particleId;
          sp = seed->sp()[1];
          std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                    << sp->z() << ") hit_id? = " << sp->particles[0].particleId;
          sp = seed->sp()[2];
          std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                    << sp->z() << ") hit_id = " << sp->particles[0].particleId
                    << " with count " << sp->particles[0].hitCount;
          std::cout << std::endl;
          printCounter2++;
        }
      }
    }
    // prints first few space points
    if (printCounter < maxNumToPrint) {
      for (size_t i = 0; i < regionVec.size(); i++) {
        const Acts::Seed<SpacePoint>* seed = &regionVec[i];
        const SpacePoint* sp = seed->sp()[0];
        std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z()
                  << ") hit_id? = " << sp->particles[0].particleId.particle();
        sp = seed->sp()[1];
        std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z()
                  << ") hit_id? = " << sp->particles[0].particleId.particle();
        sp = seed->sp()[2];
        std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z()
                  << ") hit_id = " << sp->particles[0].particleId.particle()
                  << " with count " << sp->particles[0].hitCount;
        std::cout << std::endl;
        printCounter++;
      }
    }
  }
  ACTS_INFO("Number of seeds generated: " << numSeeds)
  ACTS_INFO("Number of true seeds generated: " << numTrueSeeds)
  ACTS_INFO("Number of empty regions: " << emptySeedCount)
  ACTS_INFO("Number of clusters (hits) used is: " << clustCounter)
  ACTS_INFO("Number of clusters (hits) ignored is: " << numIgnored)

  return FW::ProcessCode::SUCCESS;
}
