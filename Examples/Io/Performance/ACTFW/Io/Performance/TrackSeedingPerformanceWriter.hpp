// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

//#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Validation/EffPlotTool.hpp"
#include "ACTFW/Validation/ResPlotTool.hpp"
#include "ACTFW/Validation/TrackSummaryPlotTool.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SpacePoint.hpp"

#include <mutex>

class TFile;
class TTree;

namespace FW {

/// Write out the residual and pull of track parameters and efficiency.
///
/// Efficiency here is the fraction of smoothed tracks compared to all tracks.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class TrackSeedingPerformanceWriter final
    : public WriterT<std::vector<std::vector<Acts::Seed<SpacePoint>>>> {
 public:
  struct Config {
    /// Input truth particles collection.
    std::string inputParticles;
    /// Input (fitted) trajectories collection.
    std::string inputSeeds;
    /// Output directory.
    std::string outputDir;
    /// Output filename.
    std::string outputFilename = "performance_track_seeding.root";
    /// Plot tool configurations.
    EffPlotTool::Config effPlotToolConfig;
  };

  /// Finds whether or not the seed contains a truth particle.
  /// Technically, space points can have multiple particles that are a part of
  /// them, so analyzeSeed finds how many particles are in common.
  /// @param seed The seed to be processed.
  /// @param particlesFoundBySeeds The set of particle barcodes already found.
  /// @param nDuplicateSeeds Number of seeds that find a particle already
  /// identified by a previous seed.
  ///
  /// Returns the number particles that are a part
  /// of all 3 spacePoints in the seed. Returning 0 means it's a fake seed.
  std::size_t analyzeSeed(const Acts::Seed<SpacePoint>* seed,
                          std::set<ActsFatras::Barcode>& particlesFoundBySeeds,
                          std::size_t& nDuplicateSeeds) const;

  void printSeed(const Acts::Seed<SpacePoint>* seed) const;

  /// Construct from configuration and log level.
  TrackSeedingPerformanceWriter(Config cfg, Acts::Logging::Level lvl);
  ~TrackSeedingPerformanceWriter() override;

  /// Finalize plots.
  ProcessCode endRun() final override;

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const std::vector<std::vector<Acts::Seed<SpacePoint>>>&
                         seedVector) final override;

  Config m_cfg;
  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Plot tool for efficiency
  EffPlotTool m_effPlotTool;
  EffPlotTool::EffPlotCache m_effPlotCache;
};

}  // namespace FW