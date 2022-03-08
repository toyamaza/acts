// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Digitization/Segmentation.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.h"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

/// @class SpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the clusters pixel detectors need further treatment. This class takes
/// the digitized clusters on a pixel detector element and provides the
/// corresponding space point.
///
template <typename spacepoint_t>
class SpacePointBuilder {
 public:
  using Measurement = Acts::BoundVariantMeasurement;
  // Constructor
  SpacePointBuilder(SpacePointBuilderConfig cfg);
  // Default constructor
  SpacePointBuilder() = default;

  void calculateSpacePoints(
      const GeometryContext& gctx, std::vector<spacepoint_t>& spacePointStorage,
      const std::vector<const Measurement>* frontMeasurements,
      const std::vector<const Measurement>* backMeasurements = nullptr) const;

 protected:
  /// @brief Getter method for the local coordinates of a cluster
  /// on its corresponding surface
  ///
  /// @param clus cluster that holds the neccesary information of the 2D hit position.
  /// @return vector of the local coordinates of the cluster on the surface
  Vector2 getLocalPos(const Measurement& meas) const;
  std::pair<Acts::Vector2, Acts::SymMatrix2> getLocalPosCov(
      const Measurement& meas) const;

  /// @brief Getter method for the global coordinates of a cluster
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param cluster cluster that holds the necessary
  /// information
  /// @return vector of the global coordinates and covariance of the cluster
  std::pair<Vector3, Vector2> globalCoords(const GeometryContext& gctx,
                                           const Measurement& meas) const;

  /// @brief Calculates the space points out of a given collection of clusters
  /// and stores the results
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param clusters vector of clusters
  /// @param spacePointStorage storage of the results
  void calculateSingleHitSpacePoints(
      const GeometryContext& gctx,
      const std::vector<const Measurement>& measurements,
      std::vector<spacepoint_t>& spacePointStorage) const;

  /// @brief Searches possible combinations of two clusters on different
  /// surfaces that may come from the same particles
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param clustersFront vector of clusters on a surface
  /// @param clustersBack vector of clusters on another surface
  /// @param clusterPairs storage of the cluster pairs
  /// @note The structure of @p clustersFront and @p clustersBack is
  /// meant to be clusters[Independent clusters on a single surface]
  void makeMeasurementPairs(
      const GeometryContext& gctx,
      const std::vector<const Measurement>& measurementsFront,
      const std::vector<const Measurement>& measurementsBack,
      std::vector<std::pair<const Measurement*, const Measurement*>>&
          measurementPairs) const;

  void calculateDoubleHitSpacePoints(
      const Acts::GeometryContext& gctx,
      const std::vector<std::pair<const Measurement*, const Measurement*>>&
          measurementPairs,
      std::vector<spacepoint_t>& spacePoints) const;

  std::pair<Acts::Vector3, Acts::Vector3> endsOfStrip(
      const Acts::GeometryContext& gctx, const Measurement& measurement) const;

  // Acts::SourceLink getSourceLink( const Measurement& meas) const;

  Acts::Vector2 globalCov(const Acts::GeometryContext& gctx,
                          const Acts::GeometryIdentifier& geoId,
                          const Acts::Vector2& localPos,
                          const Acts::SymMatrix2& localCov) const;

  double getLocVar(const Measurement& meas) const;

  Acts::Vector2 getGlobalVars(const Acts::GeometryContext& gctx,
                              const Measurement& meas_front,
                              const Measurement& meas_back,
                              const double theta) const;

  const Acts::SourceLink* getSourceLink(const Measurement meas) const;

  // configuration of the single hit space point builder
  SpacePointBuilderConfig m_config;
  // Acts::Vector3 globalPos(const Acts::GeometryContext& gctx, const
  // Measurement& meas) const;
};

}  // namespace Acts
#include "Acts/SpacePointFormation/detail/SpacePointBuilder.ipp"
