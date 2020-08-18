// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/GenericDetector/GenericDetector.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Csv/CsvOptionsReader.hpp"
#include "ACTFW/Io/Csv/CsvOptionsWriter.hpp"
#include "ACTFW/Io/Csv/CsvParticleReader.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterWriter.hpp"
#include "ACTFW/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ACTFW/Io/Performance/TrackSeedingPerformanceWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/TruthTracking/TruthSeedSelector.hpp"
#include "ACTFW/TruthTracking/TruthTrackFinder.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include <Acts/Utilities/Units.hpp>

#include <cstdlib>
#include <memory>

#include "TestSeedAlgorithm.hpp"

using namespace FW;

int main(int argc, char* argv[]) {
  GenericDetector detector;

  // setup options
  // every component should have an associated option setup function
  // that should be called here.
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  detector.addOptions(desc);
  Options::addBFieldOptions(desc);
  Options::addCsvWriterOptions(desc);

  // parse options from command line flags
  auto vm = Options::parse(desc, argc, argv);
  // an empty varaibles map indicates an error
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // setup the sequencer first w/ config derived from options
  Sequencer sequencer(FW::Options::readSequencerConfig(vm));

  // extract some common options
  auto logLevel = Options::readLogLevel(vm);
  auto inputDir = vm["input-dir"].as<std::string>();
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());

  // Setup detector geometry
  auto geometry = Geometry::build(vm, detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup the magnetic field
  auto magneticField = Options::readBField(vm);
  /*
    // Read particles (initial states) and clusters from CSV files
    auto particleReader = Options::readCsvParticleReaderConfig(vm);
    particleReader.inputStem = "particles_initial";
    particleReader.outputParticles = "particles_initial";
    // sequencer.addReader(
    //  std::make_shared<CsvParticleReader>(particleReader, logLevel));
    */

  // Read particles (initial states) and clusters from CSV files
  auto particleReader = Options::readCsvParticleReaderConfig(vm);
  particleReader.inputStem = "particles_initial";
  particleReader.outputParticles = "particles_initial";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(particleReader, logLevel));
  // Read clusters from CSV files
  auto clusterReaderCfg = Options::readCsvPlanarClusterReaderConfig(vm);
  clusterReaderCfg.trackingGeometry = trackingGeometry;
  clusterReaderCfg.outputClusters = "clusters";
  clusterReaderCfg.outputHitIds = "hit_ids";
  clusterReaderCfg.outputHitParticlesMap = "hit_particles_map";
  clusterReaderCfg.outputSimulatedHits = "hits";
  sequencer.addReader(
      std::make_shared<CsvPlanarClusterReader>(clusterReaderCfg, logLevel));

  const auto& inputParticles = particleReader.outputParticles;

  // Cuts down on the number of particles so that effeciency is calculated
  // properly
  TruthSeedSelector::Config seedSelectorCfg;
  seedSelectorCfg.inputParticles = inputParticles;
  seedSelectorCfg.inputHitParticlesMap = clusterReaderCfg.outputHitParticlesMap;
  seedSelectorCfg.outputParticles =
      "particles_final";  // not sure this is what it should be called.
                          // particles_final may already have a different
                          // meaning
  seedSelectorCfg.etaMin = -2.7;
  seedSelectorCfg.etaMax = 2.7;
  seedSelectorCfg.ptMin = 0.5;
  seedSelectorCfg.nHitsMin = 3;
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(seedSelectorCfg, logLevel));

  // add Seeding Algorithm that finds the seeds
  FW::TestSeedAlgorithm::Config testSeedCfg;
  testSeedCfg.inputHitParticlesMap = "hit_particles_map";
  testSeedCfg.inputSimulatedHits = "hits";
  testSeedCfg.inputDir = inputDir;
  testSeedCfg.outputHitIds = "hit_ids";
  testSeedCfg.inputClusters = "clusters";
  testSeedCfg.outputSeeds = "output_seeds";
  sequencer.addAlgorithm(
      std::make_shared<FW::TestSeedAlgorithm>(testSeedCfg, logLevel));

  FW::TrackSeedingPerformanceWriter::Config seedPerfCfg;
  seedPerfCfg.inputSeeds = testSeedCfg.outputSeeds;
  seedPerfCfg.inputParticles = seedSelectorCfg.outputParticles;
  seedPerfCfg.outputDir = outputDir;
  sequencer.addWriter(
      std::make_shared<TrackSeedingPerformanceWriter>(seedPerfCfg, logLevel));
  // Run all configured algorithms and return the appropriate status.
  return sequencer.run();
}
