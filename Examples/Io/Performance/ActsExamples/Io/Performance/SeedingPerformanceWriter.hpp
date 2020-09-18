// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Seeding/SimSpacePoint.hpp"
#include "ActsExamples/Validation/DuplicationPlotTool.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/FakeRatePlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <mutex>

class TFile;
class TTree;

// using namespace Acts::UnitLiterals;

namespace ActsExamples {

class SeedingPerformanceWriter final
    : public WriterT<std::vector<std::vector<Acts::Seed<SimSpacePoint>>>> {
 public:
  struct Config {
    /// Input hit to particles map
    std::string inputHitParticlesMap;
    /// Input truth particles collection.
    std::string inputParticles;
    /// Input seeds to be analyzed.
    std::string inputSeeds;
    /// Input seeds as proto tracks
    std::string inputProtoTracks;
    /// input Clusters from the event#-hits.csv file for calculating efficiency.
    std::string inputClusters;
    /// Output directory.
    std::string outputDir;
    /// Output filename.
    std::string outputFilename = "performance_track_seeding.root";

    // The quality cuts to be applied when evaluating seed finder efficiency
    /// Maximum distance from the origin in the transverse plane
    double rhoMax = std::numeric_limits<double>::max();
    /// Maximum absolute distance from the origin along z
    double absZMax = std::numeric_limits<double>::max();
    // Truth particle kinematic cuts
    double phiMin = std::numeric_limits<double>::lowest();
    double phiMax = std::numeric_limits<double>::max();
    double etaMin = std::numeric_limits<double>::lowest();
    double etaMax = std::numeric_limits<double>::max();
    double absEtaMin = std::numeric_limits<double>::lowest();
    double absEtaMax = std::numeric_limits<double>::max();
    double ptMin = 0.0;
    double ptMax = std::numeric_limits<double>::max();
    /// Keep neutral particles
    bool keepNeutral = false;
    /// Requirement on number of recorded hits
    //@TODO: implement detector-specific requirements
    size_t nHitsMin = 0;
    size_t nHitsMax = std::numeric_limits<size_t>::max();

    /// Plot tool configurations.
    EffPlotTool::Config effPlotToolConfig;
    FakeRatePlotTool::Config fakeRatePlotConfig;
  };

  /// @brief Finds all the particles that are in common to all space points in
  /// the seed.
  /// @param seed The seed to be analyzed
  std::set<ActsFatras::Barcode> identifySharedParticles(
      const Acts::Seed<SimSpacePoint>* seed) const;

  /// @brief Returns true if we expect the seed finder to be able to find this
  /// particle
  /// @param prt The particle to see whether it's findable
  /// @param particleHitsMap inverted map from hitParticlesMap, maps from
  /// particle barcodes to hits. TODO: get type information
  /// @param clusters Used to get information on the hits that a particle has
  /// made. TODO: Verify this is a valid filter
  bool prtFindable(
      const ActsFatras::Particle& prt,
      // const auto& particleHitsMap,
      const IndexMultimap<unsigned long,ActsFatras::Barcode>&  particleHitsMap,
      const ActsExamples::GeometryIdMultimap<Acts::PlanarModuleCluster>&
          clusters) const;

  /// @brief Analyzes onse seed. Finds whether or not the seed contains a truth
  /// particle. Technically, space points can have multiple particles that are a
  /// part of them, so analyzeSeed finds how many particles are in common.
  /// @param seed The seed to be processed.
  /// @param hitParticlesMap map from hits to particles
  /// @param truthCount map from particles found to how many seeds found them
  /// @param fakeCount map from particles to how many fake seeds they were a
  /// part of
  void analyzeSeed(
      const Acts::Seed<SimSpacePoint>* seed,
      const ActsExamples::IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
      std::unordered_map<ActsFatras::Barcode, std::size_t>& truthCount,
      std::unordered_map<ActsFatras::Barcode, std::size_t>& fakeCount) const;

  /// @brief Converts a seed to a proto track of 3 hits
  ActsExamples::ProtoTrack seedToProtoTrack(
      const Acts::Seed<SimSpacePoint>* seed) const;

  /// @brief Writes the fake rate and efficiency plots
  /// TODO: Add the rests of the fake rate plots.
  void writePlots(
      const std::vector<std::vector<Acts::Seed<SimSpacePoint>>>& seedVector,
      const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
      const SimParticleContainer& particles,
      const ActsExamples::GeometryIdMultimap<Acts::PlanarModuleCluster>&
          clusters);

  /// Construct from configuration and log level.
  SeedingPerformanceWriter(Config cfg, Acts::Logging::Level lvl);
  ~SeedingPerformanceWriter() override;

  /// Finalize plots.
  ProcessCode endRun() final override;

 private:
  /// @brief Calls TTrees writing function, and the writePlots function
  /// These create two files, one that has premade plots, and another that can
  /// be browsed.
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const std::vector<std::vector<Acts::Seed<SimSpacePoint>>>&
                         seedVector) final override;
  /// @brief Used for writing the TTrees
  struct Impl;
  std::unique_ptr<Impl> m_impl;

  Config m_cfg;
  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Plot tool for fakerate
  FakeRatePlotTool m_fakeRatePlotTool;
  FakeRatePlotTool::FakeRatePlotCache m_fakeRatePlotCache;
  /// Plot tool for efficiency
  EffPlotTool m_effPlotTool;
  EffPlotTool::EffPlotCache m_effPlotCache;
};

}  // namespace ActsExamples
