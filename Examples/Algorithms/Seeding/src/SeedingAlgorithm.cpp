// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Seeding/SeedingAlgorithm.hpp"

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/EventData/SimHit.hpp"
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

#include "ACTFW/Seeding/SpacePointFromHit.hpp"
#include "ACTFW/Seeding/GenericDetectorCuts.hpp"

#include <iostream>
#include <stdexcept>

FW::SeedingAlgorithm::SeedingAlgorithm(
    FW::SeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : FW::BareAlgorithm("SeedingAlgorithm", lvl), m_cfg(std::move(cfg)) {
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

FW::ProcessCode FW::SeedingAlgorithm::execute(const AlgorithmContext& ctx) const {

  // Prepare the input and output collections
  const auto& hits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);

  Acts::SeedfinderConfig<SpacePointFromHit> config;
  // silicon detector max
  // config.rMax = 160.;
  // config.rMax = 300.;
  // config.deltaRMin = 5.;
  // config.deltaRMax = 160.;
  // config.collisionRegionMin = -250.;
  // config.collisionRegionMax = 250.;
  // config.zMin = -2800.;
  // config.zMax = 2800.;
  // config.maxSeedsPerSpM = 5;
  // config.cotThetaMax = 7.40627;  // 2.7 eta
  // config.sigmaScattering = 1.00000;
  // config.minPt = 500.;
  // config.bFieldInZ = 0.00199724;
  // config.beamPos = {-.5, -.5};
  // config.impactMax = 10.;

  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePointFromHit>>(
									      Acts::BinFinder<SpacePointFromHit>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SpacePointFromHit>>(
									   Acts::BinFinder<SpacePointFromHit>());
  Acts::SeedFilterConfig sfconf;
  Acts::GenericDetectorCuts<SpacePointFromHit> atlasCuts = Acts::GenericDetectorCuts<SpacePointFromHit>();
  config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePointFromHit>>(
									    Acts::SeedFilter<SpacePointFromHit>(sfconf, &atlasCuts));
  Acts::Seedfinder<SpacePointFromHit> seedFinder(config);


  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SpacePointFromHit& sp, float, float, float) -> Acts::Vector2D {
    return {sp.varianceR, sp.varianceZ};
  };


  std::vector<const SpacePointFromHit*> spVec;
  for (const auto& hit : hits) {
    const auto hitPos = hit.position();
    int layer = 1; // dummy
    float varR = 0.01;
    float varZ = 0.5;
    float hitPosX = hitPos.x();
    float hitPosY = hitPos.y();
    float hitPosZ = hitPos.z();
    float r = sqrt(hitPosX*hitPosX + hitPosY*hitPosY);
    SpacePointFromHit* sp = new SpacePointFromHit{
      hitPosX,hitPosY,hitPosZ, r, layer, varR, varZ
    };
    spVec.push_back(sp);
  }


  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<Acts::SpacePointGrid<SpacePointFromHit>> grid =
    Acts::SpacePointGridCreator::createGrid<SpacePointFromHit>(gridConf);
  auto spGroup = Acts::BinnedSPGroup<SpacePointFromHit>(spVec.begin(), spVec.end(), ct,
							bottomBinFinder, topBinFinder,
							std::move(grid), config);


  std::vector<std::vector<Acts::Seed<SpacePointFromHit>>> seedVector;
  // auto start = std::chrono::system_clock::now();
  auto groupIt = spGroup.begin();
  auto endOfGroups = spGroup.end();
  int cnt = 0;
  for (; !(groupIt == endOfGroups); ++groupIt) {
    // std::cout << "group" << cnt << std::endl;
    cnt++;

    seedVector.push_back(seedFinder.createSeedsForGroup(
							groupIt.bottom(), groupIt.middle(), groupIt.top()));
  }
    
  int numSeeds = 0;
  for (auto& outVec : seedVector) {
    numSeeds += outVec.size();
  }

  std::cout << spVec.size() << " hits, " << seedVector.size() << " regions, " << numSeeds << " seeds" << std::endl;

  for (auto& regionVec : seedVector) {
    for (size_t i = 0; i < regionVec.size(); i++) {
      const Acts::Seed<SpacePointFromHit>* seed = &regionVec[i];
      const SpacePointFromHit* sp = seed->sp()[0];
      std::cout << " (" << sp->x() << ", " << sp->y() << ", " << sp->z()
		<< ") ";
      sp = seed->sp()[1];
      std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
		<< sp->z() << ") ";
      sp = seed->sp()[2];
      std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
		<< sp->z() << ") ";
      std::cout << std::endl;
    }
  }
  std::cout << std::endl;

  return FW::ProcessCode::SUCCESS;
}
