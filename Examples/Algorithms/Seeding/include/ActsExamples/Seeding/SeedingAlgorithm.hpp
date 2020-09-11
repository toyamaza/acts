// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples//Framework/BareAlgorithm.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/Seeding/SimSpacePoint.hpp"


#include <set>
#include <memory>
#include <string>
#include <unordered_map>

namespace Acts {
// class DigitizationModule;
class IdentifiedDetectorElement;
class PlanarModuleStepper;
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace ActsExamples {

/// Create planar clusters from simulation hits.
class SeedingAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input collection of simulated hits.
    std::string inputSimulatedHits;
    /// Output collection of clusters.
    std::string outputSeeds;
    /// Output prototracks
    std::string outputProtoTracks;
    // input Clusters from the event#-hits.csv file.
    std::string inputClusters;
    // input particles for creating proto seeds
    std::string inputParticles;
    // input dir containing hits and truth information
    std::string inputDir;
    // used to get truth information into seeds about what particles are in what
    // space point.
    std::string inputHitParticlesMap;
    /// Which simulated (truth) hits collection to use. Not used currently.

    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

  };

  /// Construct the digitization algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Technically, space points can have multiple particles that are a part of
  /// them, so analyzeSeed finds how many particles are in common.
  /// @param seed The seed to be processed.
  /// @param particlesFoundBySeeds The set of particle barcodes already found.
  /// @param nDuplicateSeeds Number of seeds that find a particle already
  /// identified by a previous seed.
  ///
  /// Returns the number particles that are a part
  /// of all 3 spacePoints in the seed. Returning 0 means it's a fake seed.
  std::size_t analyzeSeed(const Acts::Seed<ActsExamples::SimSpacePoint>* seed,
                          std::set<ActsFatras::Barcode>& particlesFoundBySeeds,
                          std::size_t& nDuplicateSeeds) const;

  /// @param cluster The hit with local hit information
  /// Returns a space point with a particle barcode stored in .particles for
  /// each particle that made this space point.
  ActsExamples::SimSpacePoint* readSP(std::size_t hit_id, const Acts::GeometryID geoId,
                     const Acts::PlanarModuleCluster& cluster,
                     const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
                     const AlgorithmContext& ctx) const;


  /// @brief Converts a seed to a proto track of 3 hits
  ActsExamples::ProtoTrack seedToProtoTrack(const Acts::Seed<ActsExamples::SimSpacePoint>* seed) const;
  
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  // struct Digitizable {
  //   const Acts::Surface* surface = nullptr;
  //   const Acts::IdentifiedDetectorElement* detectorElement = nullptr;
  //   const Acts::DigitizationModule* digitizer = nullptr;
  // };

  Config m_cfg;
  /// Lookup container for all digitizable surfaces
  std::unordered_map<Acts::GeometryID, const Acts::Surface*> m_surfaces;
  // std::unordered_map<Acts::GeometryID, Digitizable> m_digitizables;
};

}  // namespace ActsExamples
