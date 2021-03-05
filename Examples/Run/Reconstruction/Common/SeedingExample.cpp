// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Digitization/HitSmearing.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedingOptions.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Io/Json/JsonSeedingConfig.hpp"

#include <memory>

#include <boost/program_options.hpp>

#include "RecInput.hpp"

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int runSeedingExample(int argc, char* argv[],
                      std::shared_ptr<ActsExamples::IBaseDetector> detector,
                      std::vector<Acts::GeometryIdentifier> layersForSeeding) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  Options::addInputOptions(desc);
  Options::addMagneticFieldOptions(desc);
  // Add specific options for this geometry
  detector->addOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }
  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Now read the standard options
  auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());

  // The geometry, material and decoration
  auto geometry = Geometry::build(vm, *detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  auto rnd =
    std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vm));

  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // Setup the magnetic field
  Options::setupMagneticFieldServices(vm, sequencer);
  auto magneticField = Options::readMagneticField(vm);

  // Read the sim hits
  auto simHitReaderCfg = setupSimHitReading(vm, sequencer);
  // Read the particles
  auto particleReader = setupParticleReading(vm, sequencer);

  // Run the sim hits smearing
  auto hitSmearingCfg = setupSimHitSmearing(vm, sequencer, rnd, tGeometry,
                                            simHitReaderCfg.outputSimHits);
  // Read json config file
  auto cfile = vm["seed-config-file"].as<std::string>();
  
  if (not cfile.empty()) {

    auto in = std::ifstream(cfile, std::ifstream::in | std::ifstream::binary);
    if (in.good()) {
      // Get the json file for the configuration
      nlohmann::json djson;
      in >> djson;
      in.close();

    }
  }
  // Run the particle selection
  // The pre-selection will select truth particles satisfying provided criteria
  // from all particles read in by particle reader for further processing. It
  // has no impact on the truth hits read-in by the cluster reader.
  TruthSeedSelector::Config particleSelectorCfg;
  particleSelectorCfg.inputParticles = particleReader.outputParticles;
  particleSelectorCfg.inputMeasurementParticlesMap =
    hitSmearingCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.ptMin = 500_MeV;
  particleSelectorCfg.etaMax = 2.5;
  particleSelectorCfg.etaMin = -2.5;
  particleSelectorCfg.nHitsMin = 9;
  sequencer.addAlgorithm(
			 std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The selected particles
  const auto& inputParticles = particleSelectorCfg.outputParticles;

  // Create space points
  SpacePointMaker::Config spCfg;
  spCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  spCfg.inputMeasurements = hitSmearingCfg.outputMeasurements;
  spCfg.outputSpacePoints = "spacepoints";
  spCfg.trackingGeometry = tGeometry;
  spCfg.geometrySelection = layersForSeeding;
  sequencer.addAlgorithm(std::make_shared<SpacePointMaker>(spCfg, logLevel));

  // Seeding algorithm
  SeedingAlgorithm::Config seedingCfg;
  seedingCfg.inputSpacePoints = {
    spCfg.outputSpacePoints,
  };
  seedingCfg.outputSeeds = "seeds";
  seedingCfg.outputProtoTracks = "prototracks";
  seedingCfg.rMax = 200.;
  seedingCfg.deltaRMax = 60.;
  seedingCfg.collisionRegionMin = -250;
  seedingCfg.collisionRegionMax = 250.;
  seedingCfg.zMin = -2000.;
  seedingCfg.zMax = 2000.;
  seedingCfg.maxSeedsPerSpM = 1;
  seedingCfg.cotThetaMax = 7.40627;  // 2.7 eta
  seedingCfg.sigmaScattering = 50;
  seedingCfg.radLengthPerSeed = 0.1;
  seedingCfg.minPt = 500.;
  seedingCfg.bFieldInZ = 0.00199724;
  seedingCfg.beamPosX = 0;
  seedingCfg.beamPosY = 0;
  seedingCfg.impactMax = 3.;
  sequencer.addAlgorithm(
			 std::make_shared<SeedingAlgorithm>(seedingCfg, logLevel));

  // Algorithm estimating track parameter from seed
  TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
  paramsEstimationCfg.inputSeeds = seedingCfg.outputSeeds;
  paramsEstimationCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
  paramsEstimationCfg.outputTrackParametersSeedMap = "estimatedparams_seed_map";
  paramsEstimationCfg.trackingGeometry = tGeometry;
  paramsEstimationCfg.magneticField = magneticField;
  sequencer.addAlgorithm(std::make_shared<TrackParamsEstimationAlgorithm>(
									  paramsEstimationCfg, logLevel));

  // Seeding performance Writers
  TrackFinderPerformanceWriter::Config tfPerfCfg;
  tfPerfCfg.inputProtoTracks = seedingCfg.outputProtoTracks;
  tfPerfCfg.inputParticles = inputParticles;
  tfPerfCfg.inputMeasurementParticlesMap =
    hitSmearingCfg.outputMeasurementParticlesMap;
  tfPerfCfg.outputDir = outputDir;
  tfPerfCfg.outputFilename = "performance_seeding_trees.root";
  sequencer.addWriter(
		      std::make_shared<TrackFinderPerformanceWriter>(tfPerfCfg, logLevel));

  SeedingPerformanceWriter::Config seedPerfCfg;
  seedPerfCfg.inputSeeds = seedingCfg.outputSeeds;
  seedPerfCfg.inputParticles = inputParticles;
  seedPerfCfg.inputMeasurementParticlesMap =
    hitSmearingCfg.outputMeasurementParticlesMap;
  seedPerfCfg.outputDir = outputDir;
  seedPerfCfg.outputFilename = "performance_seeding_hists.root";
  sequencer.addWriter(
		      std::make_shared<SeedingPerformanceWriter>(seedPerfCfg, logLevel));

  // The track parameters estimation writer
  RootTrackParameterWriter::Config trackParamsWriterCfg;
  trackParamsWriterCfg.inputSeeds = seedingCfg.outputSeeds;
  trackParamsWriterCfg.inputTrackParameters =
    paramsEstimationCfg.outputTrackParameters;
  trackParamsWriterCfg.inputTrackParametersSeedMap =
    paramsEstimationCfg.outputTrackParametersSeedMap;
  trackParamsWriterCfg.inputParticles = particleReader.outputParticles;
  trackParamsWriterCfg.inputSimHits = simHitReaderCfg.outputSimHits;
  trackParamsWriterCfg.inputMeasurementParticlesMap =
    hitSmearingCfg.outputMeasurementParticlesMap;
  trackParamsWriterCfg.inputMeasurementSimHitsMap =
    hitSmearingCfg.outputMeasurementSimHitsMap;
  trackParamsWriterCfg.outputDir = outputDir;
  trackParamsWriterCfg.outputFilename = "estimatedparams.root";
  trackParamsWriterCfg.outputTreename = "estimatedparams";
  sequencer.addWriter(std::make_shared<RootTrackParameterWriter>(
								 trackParamsWriterCfg, logLevel));

  return sequencer.run();
}
