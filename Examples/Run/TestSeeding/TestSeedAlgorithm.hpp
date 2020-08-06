// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "ACTFW/EventData/IndexContainers.hpp"
#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SpacePoint.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

namespace FW {

/// A simple empty algorithm
class TestSeedAlgorithm : public FW::BareAlgorithm {
 public:
  struct Config {
    // currently not using outputHitIds.
    std::string outputHitIds;
    // input Clusters from the event#-hits.csv file.
    std::string inputClusters;
    // inputDir not currently used.
    std::string inputDir;
    // used to get truth information into seeds about what particles are in what
    // space point.
    std::string inputHitParticlesMap;
    /// Which simulated (truth) hits collection to use. Not used currently.
    std::string inputSimulatedHits;
  };

  TestSeedAlgorithm(const Config& cfg, Acts::Logging::Level level);

  std::size_t seedNumParticles(const Acts::Seed<SpacePoint>* seed) const;

  SpacePoint* readSP(std::size_t hit_id, const Acts::GeometryID geoId,
                     const Acts::PlanarModuleCluster& cluster,
                     const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
                     const AlgorithmContext& ctx) const;
  /// The framework execut mehtod
  /// @param ctx The Algorithm context for multithreading
  FW::ProcessCode execute(
      const AlgorithmContext& ctx /*,
      const GeometryIdMultimap<Acts::PlanarModuleCluster>& clusters*/) const final override;

 private:
  Config m_cfg;
};

}  // namespace FW
