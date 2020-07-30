// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
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
}

FW::ProcessCode FW::SeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {

  // Prepare the input and output collections
  const auto& hits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);
  Acts::SeedfinderConfig<SimHit> config;
  // silicon detector max
  config.rMax = 160.;
  config.deltaRMin = 5.;
  config.deltaRMax = 160.;
  config.collisionRegionMin = -250.;
  config.collisionRegionMax = 250.;
  config.zMin = -2800.;
  config.zMax = 2800.;
  config.maxSeedsPerSpM = 5;
  // 2.7 eta
  config.cotThetaMax = 7.40627;
  config.sigmaScattering = 1.00000;

  config.minPt = 500.;
  config.bFieldInZ = 0.00199724;

  config.beamPos = {-.5, -.5};
  config.impactMax = 10.;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SimHit>>(
      Acts::BinFinder<SimHit>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SimHit>>(
      Acts::BinFinder<SimHit>());
  Acts::SeedFilterConfig sfconf;
  // Acts::ATLASCuts<SimHit> atlasCuts = Acts::ATLASCuts<SimHit>();
  // config.seedFilter = std::make_unique<Acts::SeedFilter<SimHit>>(
  //     Acts::SeedFilter<SimHit>(sfconf, &atlasCuts));
  Acts::Seedfinder<SimHit> a(config);

  covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SimHit& sp, float, float, float) -> Acts::Vector2D {
return {0., 0.}; // Fix me (cov R, cov Z)
  };

  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;
  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<Acts::SpacePointGrid<SimHit>> grid =
      Acts::SpacePointGridCreator::createGrid<SpacePoint>(gridConf);
  auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(), spVec.end(), ct,
                                                 bottomBinFinder, topBinFinder,
                                                 std::move(grid), config);

  // std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector;
  // auto start = std::chrono::system_clock::now();
  // auto groupIt = spGroup.begin();
  // auto endOfGroups = spGroup.end();
  // for (; !(groupIt == endOfGroups); ++groupIt) {
  //   seedVector.push_back(a.createSeedsForGroup(
  //       groupIt.bottom(), groupIt.middle(), groupIt.top()));
  // }
  // auto end = std::chrono::system_clock::now();
  // std::chrono::duration<double> elapsed_seconds = end - start;
  // std::cout << "time to create seeds: " << elapsed_seconds.count() << std::endl;
  // std::cout << "Number of regions: " << seedVector.size() << std::endl;
  // int numSeeds = 0;
  // for (auto& outVec : seedVector) {
  //   numSeeds += outVec.size();
  // }
  // std::cout << "Number of seeds generated: " << numSeeds << std::endl;

// for (auto&& [moduleGeoId, moduleHits] : groupByModule(hits)) {
  //   // can only digitize hits on digitizable surfaces
  //   const auto it = m_digitizables.find(moduleGeoId);
  //   if (it == m_digitizables.end()) {
  //     continue;
  //   }

  //   const auto& dg = it->second;
  //   // local intersection / direction
  //   const auto invTransfrom = dg.surface->transform(ctx.geoContext).inverse();

  //   // use iterators manually so we can retrieve the hit index in the container
  //   for (auto ih = moduleHits.begin(); ih != moduleHits.end(); ++ih) {
  //     const auto& hit = *ih;
  //     const auto idx = hits.index_of(ih);

  //     Acts::Vector2D localIntersect = (invTransfrom * hit.position()).head<2>();
  //     Acts::Vector3D localDirection =
  //         invTransfrom.linear() * hit.unitDirection();

  //     // compute digitization steps
  //     const auto thickness = dg.detectorElement->thickness();
  //     const auto lorentzAngle = dg.digitizer->lorentzAngle();
  //     auto lorentzShift = thickness * std::tan(lorentzAngle);
  //     lorentzShift *= -(dg.digitizer->readoutDirection());
  //     // now calculate the steps through the silicon
  //     std::vector<Acts::DigitizationStep> dSteps =
  //         m_cfg.planarModuleStepper->cellSteps(ctx.geoContext, *dg.digitizer,
  //                                              localIntersect, localDirection);
  //     // everything under threshold or edge effects
  //     if (!dSteps.size()) {
  //       ACTS_VERBOSE("No steps returned from stepper.");
  //       continue;
  //     }

  //     // lets create a cluster - centroid method
  //     double localX = 0.;
  //     double localY = 0.;
  //     double totalPath = 0.;
  //     // the cells to be used
  //     std::vector<Acts::DigitizationCell> usedCells;
  //     usedCells.reserve(dSteps.size());
  //     // loop over the steps
  //     for (auto dStep : dSteps) {
  //       // @todo implement smearing
  //       localX += dStep.stepLength * dStep.stepCellCenter.x();
  //       localY += dStep.stepLength * dStep.stepCellCenter.y();
  //       totalPath += dStep.stepLength;
  //       usedCells.push_back(Acts::DigitizationCell(dStep.stepCell.channel0,
  //                                                  dStep.stepCell.channel1,
  //                                                  dStep.stepLength));
  //     }
  //     // divide by the total path
  //     localX /= totalPath;
  //     localX += lorentzShift;
  //     localY /= totalPath;

  //     // get the segmentation & find the corresponding cell id
  //     const Acts::Segmentation& segmentation = dg.digitizer->segmentation();
  //     auto binUtility = segmentation.binUtility();
  //     Acts::Vector2D localPosition(localX, localY);
  //     // @todo remove unneccesary conversion
  //     // size_t bin0 = binUtility.bin(localPosition, 0);
  //     // size_t bin1 = binUtility.bin(localPosition, 1);
  //     // size_t binSerialized = binUtility.serialize({{bin0, bin1, 0}});

  //     // the covariance is currently set to 0.
  //     Acts::ActsSymMatrixD<3> cov;
  //     cov << 0.05, 0., 0., 0., 0.05, 0., 0., 0.,
  //         900. * Acts::UnitConstants::ps * Acts::UnitConstants::ps;

  //     // create the planar cluster
  //     Acts::PlanarModuleCluster pCluster(
  //         dg.surface->getSharedPtr(), Identifier(identifier_type(idx), {idx}),
  //         std::move(cov), localX, localY, hit.time(), std::move(usedCells));

  //     // insert into the cluster container. since the input data is already
  //     // sorted by geoId, we should always be able to add at the end.
  //     clusters.emplace_hint(clusters.end(), hit.geometryId(),
  //                           std::move(pCluster));
  //   }
  // }

  // ACTS_DEBUG("digitized " << hits.size() << " hits into " << clusters.size()
  //                         << " clusters");

  // // write the clusters to the EventStore
  // ctx.eventStore.add(m_cfg.outputClusters, std::move(clusters));
  return FW::ProcessCode::SUCCESS;
}
