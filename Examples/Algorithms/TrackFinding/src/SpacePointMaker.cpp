// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.h"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

ActsExamples::SpacePointMaker::SpacePointMaker(Config cfg,
                                               Acts::Logging::Level lvl)
    : BareAlgorithm("SpacePointMaker", lvl), m_cfg(std::move(cfg)) {
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source link input collection");
  }
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement input collection");
  }
  if (m_cfg.outputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point output collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (m_cfg.geometrySelection.empty()) {
    throw std::invalid_argument("Missing space point maker geometry selection");
  }
  // ensure geometry selection contains only valid inputs
  for (const auto& geoId : m_cfg.geometrySelection) {
    if ((geoId.approach() != 0u) or (geoId.boundary() != 0u) or
        (geoId.sensitive() != 0u)) {
      throw std::invalid_argument(
          "Invalid geometry selection: only volume and layer are allowed to be "
          "set");
    }
  }
  // remove geometry selection duplicates
  //
  // the geometry selections must be mutually exclusive, i.e. if we have a
  // selection that contains both a volume and a layer within that same volume,
  // we would create the space points for the layer twice.
  auto isDuplicate = [](Acts::GeometryIdentifier ref,
                        Acts::GeometryIdentifier cmp) {
    // code assumes ref < cmp and that only volume and layer can be non-zero
    // root node always contains everything
    if (ref.volume() == 0) {
      return true;
    }
    // unequal volumes always means separate hierarchies
    if (ref.volume() != cmp.volume()) {
      return false;
    }
    // within the same volume hierarchy only consider layers
    return (ref.layer() == cmp.layer());
  };
  auto geoSelBeg = m_cfg.geometrySelection.begin();
  auto geoSelEnd = m_cfg.geometrySelection.end();
  // sort geometry selection so the unique filtering works
  std::sort(geoSelBeg, geoSelEnd);
  auto geoSelLastUnique = std::unique(geoSelBeg, geoSelEnd, isDuplicate);
  if (geoSelLastUnique != geoSelEnd) {
    ACTS_WARNING("Removed " << std::distance(geoSelLastUnique, geoSelEnd)
                            << " geometry selection duplicates");
    m_cfg.geometrySelection.erase(geoSelLastUnique, geoSelEnd);
  }
  ACTS_INFO("Space point geometry selection:");
  for (const auto& geoId : m_cfg.geometrySelection) {
    ACTS_INFO("  " << geoId);
  }
  auto spBuilderConfig = Acts::SpacePointBuilderConfig();
  spBuilderConfig.trackingGeometry = m_cfg.trackingGeometry;
  m_spacePointBuilder = Acts::SpacePointBuilder<SimSpacePoint>(
      spBuilderConfig, Acts::getDefaultLogger("SpacePointBuilder", lvl));
}

ActsExamples::ProcessCode ActsExamples::SpacePointMaker::execute(
    const AlgorithmContext& ctx) const {
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);
  const auto& measurements =
      ctx.eventStore.get<MeasurementContainer>(m_cfg.inputMeasurements);

  std::vector<const Measurement*> selectedMeasurements;
  // TODO Support strip measurements

  for (Acts::GeometryIdentifier geoId : m_cfg.geometrySelection) {
    // select volume/layer depending on what is set in the geometry id
    auto range = selectLowestNonZeroGeometryObject(sourceLinks, geoId);
    // groupByModule only works with geometry containers, not with an
    // arbitrary range. do the equivalent grouping manually
    auto groupedByModule = makeGroupBy(range, detail::GeometryIdGetter());

    for (auto [moduleGeoId, moduleSourceLinks] : groupedByModule) {
      for (auto& sourceLink : moduleSourceLinks) {
        const auto& meas = measurements[sourceLink.get().index()];
        selectedMeasurements.emplace_back(&meas);
      }
    }
  }
  SimSpacePointContainer spacePoints;
  m_spacePointBuilder.calculateSpacePoints(
      ctx.geoContext, std::back_inserter(spacePoints), &selectedMeasurements);
  spacePoints.shrink_to_fit();

  ACTS_DEBUG("Created " << spacePoints.size() << " space points");
  ctx.eventStore.add(m_cfg.outputSpacePoints, std::move(spacePoints));

  return ActsExamples::ProcessCode::SUCCESS;
}
