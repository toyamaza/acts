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
#include "Acts/Seeding/Seed.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/Seeding/SimSpacePoint.hpp"


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
    
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

  };

  /// Construct the digitization algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);

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
