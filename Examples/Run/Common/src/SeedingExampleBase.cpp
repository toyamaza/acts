// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Detector/IBaseDetector.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Csv/CsvOptionsReader.hpp"
#include "ACTFW/Io/Csv/CsvParticleReader.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ACTFW/Seeding/SeedingAlgorithm.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include <Acts/Geometry/TrackingGeometry.hpp>

#include <memory>

#include <boost/program_options.hpp>

int seedingExample(int argc, char* argv[], FW::IBaseDetector& detector) {
  // Setup and parse options

  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addOutputOptions(desc);
  FW::Options::addInputOptions(desc);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }
  FW::Sequencer sequencer(FW::Options::readSequencerConfig(vm));

  // Now read the standard options
  auto logLevel = FW::Options::readLogLevel(vm);

  // The geometry, material and decoration
  auto geometry = FW::Geometry::build(vm, detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // Read particles (initial states) and clusters from CSV files
  auto particleReader = FW::Options::readCsvParticleReaderConfig(vm);
  particleReader.inputStem = "particles_initial";
  particleReader.outputParticles = "particles_initial";
  sequencer.addReader(std::make_shared<FW::CsvParticleReader>(particleReader, logLevel));

  // Read clusters from CSV files
  auto clusterReaderCfg = FW::Options::readCsvPlanarClusterReaderConfig(vm);
  clusterReaderCfg.trackingGeometry = tGeometry;
  clusterReaderCfg.outputClusters = "clusters";
  clusterReaderCfg.outputHitIds = "hit_ids";
  clusterReaderCfg.outputHitParticlesMap = "hit_particles_map";
  clusterReaderCfg.outputSimulatedHits = "hits";
  sequencer.addReader(std::make_shared<FW::CsvPlanarClusterReader>(clusterReaderCfg, logLevel));

  // Seeding algorithm
  FW::SeedingAlgorithm::Config seeding;
  seeding.inputSimulatedHits = clusterReaderCfg.outputSimulatedHits;
  seeding.outputSeeds =  "seeds";
  seeding.trackingGeometry = tGeometry;
  sequencer.addAlgorithm( std::make_shared<FW::SeedingAlgorithm>(seeding, logLevel));

 // Performance Writer
 // ... will be added later

  return sequencer.run();
}
