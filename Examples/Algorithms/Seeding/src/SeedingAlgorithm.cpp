// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Seeding/SeedingAlgorithm.hpp"

#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

#include "ActsExamples/Seeding/SimSpacePoint.hpp"
#include "ActsExamples/Seeding/GenericDetectorCuts.hpp"
#include "ActsExamples/Seeding/SeedContainer.hpp"

#include <iostream>
#include <stdexcept>

ActsExamples::SeedingAlgorithm::SeedingAlgorithm(
    ActsExamples::SeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SeedingAlgorithm", lvl), m_cfg(std::move(cfg)) {
  if (m_cfg.inputSimulatedHits.empty()) {
    throw std::invalid_argument("Missing input hits collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing output seeds collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  // fill the surface map to allow lookup by geometry id only
  m_cfg.trackingGeometry->visitSurfaces([this](const Acts::Surface* surface) {
    // for now we just require a valid surface
    if (not surface) {
      return;
    }
    this->m_surfaces.insert_or_assign(surface->geoID(), surface);
  });

}

ActsExamples::ProcessCode ActsExamples::SeedingAlgorithm::execute(const AlgorithmContext& ctx) const {

  // Prepare the input and output collections
  const auto& hits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);

  Acts::SeedfinderConfig<SimSpacePoint> config;
  // silicon detector max
  config.rMax = 200.;
  config.deltaRMin = 5.;
  config.deltaRMax = 160.;
  config.collisionRegionMin = -250;
  config.collisionRegionMax = 250.;
  config.zMin = -2000.;
  config.zMax = 2000.;
  config.maxSeedsPerSpM = 5;
  config.cotThetaMax = 7.40627;  // 2.7 eta
  config.sigmaScattering = 1.00000;
  config.minPt = 500.;
  config.bFieldInZ = 0.00199724;
  config.beamPos = {-.5, -.5};
  config.impactMax = 10.;

  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
									      Acts::BinFinder<SimSpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SimSpacePoint>>(
									   Acts::BinFinder<SimSpacePoint>());
  Acts::SeedFilterConfig sfconf;
  Acts::GenericDetectorCuts<SimSpacePoint> detectorCuts = Acts::GenericDetectorCuts<SimSpacePoint>();
  config.seedFilter = std::make_unique<Acts::SeedFilter<SimSpacePoint>>(
									    Acts::SeedFilter<SimSpacePoint>(sfconf, &detectorCuts));
  Acts::Seedfinder<SimSpacePoint> seedFinder(config);


  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SimSpacePoint& sp, float, float, float) -> Acts::Vector2D {
    return {sp.varianceR, sp.varianceZ};
  };


  std::vector<const SimSpacePoint*> spVec;

  for (const auto& hit : hits) {
    const auto hitPos = hit.position();
    int layer = 1; // dummy
    float varR = 0.01;
    float varZ = 0.5;
    float hitPosX = hitPos.x();
    float hitPosY = hitPos.y();
    float hitPosZ = hitPos.z();
    float r = sqrt(hitPosX*hitPosX + hitPosY*hitPosY);
    SimSpacePoint* sp = new SimSpacePoint{
      hitPosX,hitPosY,hitPosZ, r, layer, varR, varZ
    };

    spVec.push_back(std::move(sp));
  }

  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<Acts::SpacePointGrid<SimSpacePoint>> grid =
    Acts::SpacePointGridCreator::createGrid<SimSpacePoint>(gridConf);
  auto spGroup = Acts::BinnedSPGroup<SimSpacePoint>(spVec.begin(), spVec.end(), ct,
							bottomBinFinder, topBinFinder,
							std::move(grid), config);

  std::vector<std::vector<Acts::Seed<SimSpacePoint>>> seedVector;
  auto groupIt = spGroup.begin();
  auto endOfGroups = spGroup.end();
  for (; !(groupIt == endOfGroups); ++groupIt) {
    seedVector.push_back(seedFinder.createSeedsForGroup(groupIt.bottom(), groupIt.middle(), groupIt.top()));
  }

  SeedContainer seeds;
  int numSeeds = 0;
  for (auto& outVec : seedVector) {
    numSeeds += outVec.size();
    for (size_t i = 0; i < outVec.size(); i++) {
      const Acts::Seed<SimSpacePoint>* seed = &outVec[i];
      seeds.emplace_back(*seed);
  }
  }


  ACTS_DEBUG(spVec.size() << " hits, " << seedVector.size() << " regions, " << numSeeds << " seeds" );

  ctx.eventStore.add(m_cfg.outputSeeds, std::move(seeds) );



  return ActsExamples::ProcessCode::SUCCESS;
}
