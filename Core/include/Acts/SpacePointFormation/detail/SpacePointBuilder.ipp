// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {
namespace detail {

/// @brief Storage container for variables related to the calculation of space
/// points
struct SpacePointParameters {
  /// Vector pointing from bottom to top end of first SDE
  Vector3 q;
  /// Vector pointing from bottom to top end of second SDE
  Vector3 r;
  /// Twice the vector pointing from vertex to to midpoint of first SDE
  Vector3 s;
  /// Twice the vector pointing from vertex to to midpoint of second SDE
  Vector3 t;
  /// Cross product between SpacePointParameters::q and
  /// SpacePointParameters::s
  Vector3 qs;
  /// Cross product between SpacePointParameters::r and
  /// SpacePointParameters::t
  Vector3 rt;
  /// Magnitude of SpacePointParameters::q
  double qmag = 0.;
  /// Parameter that determines the hit position on the first SDE
  double m = 0.;
  /// Parameter that determines the hit position on the second SDE
  double n = 0.;
  /// Regular limit of the absolut values of SpacePointParameters::m and
  /// SpacePointParameters::n
  double limit = 1.;
  /// Limit of SpacePointParameters::m and SpacePointParameters::n in case of
  /// variable vertex
  double limitExtended = 0.;
};

/// @brief Calculates (Delta theta)^2 + (Delta phi)^2 between two clusters
///
/// @param [in] pos1 position of the first cluster
/// @param [in] pos2 position the second cluster
/// @param [in] maxDistance Maximum distance between two clusters
/// @param [in] maxAngleTheta2 Maximum squared theta angle between two
/// clusters
/// @param [in] maxAnglePhi2 Maximum squared phi angle between two clusters
///
/// @return The squared sum within configuration parameters, otherwise -1
inline double differenceOfMeasurementsChecked(const Acts::Vector3& pos1,
                                              const Acts::Vector3& pos2,
                                              const Acts::Vector3& posVertex,
                                              const double maxDistance,
                                              const double maxAngleTheta2,
                                              const double maxAnglePhi2) {
  std::cout << "differenceofmeasurementschecked" << std::endl;
  // Check if clusters are close enough to each other
  if ((pos1 - pos2).norm() > maxDistance) {
    return -1.;
  }

  // Calculate the angles of the vectors
  double phi1, theta1, phi2, theta2;
  phi1 = Acts::VectorHelpers::phi(pos1 - posVertex);
  theta1 = Acts::VectorHelpers::theta(pos1 - posVertex);
  phi2 = Acts::VectorHelpers::phi(pos2 - posVertex);
  theta2 = Acts::VectorHelpers::theta(pos2 - posVertex);

  // Calculate the squared difference between the theta angles
  double diffTheta2 = (theta1 - theta2) * (theta1 - theta2);
  if (diffTheta2 > maxAngleTheta2) {
    return -1.;
  }
  // Calculate the squared difference between the phi angles
  double diffPhi2 = (phi1 - phi2) * (phi1 - phi2);
  if (diffPhi2 > maxAnglePhi2) {
    return -1.;
  }
  // Return the squared distance between both vector
  return diffTheta2 + diffPhi2;
}

/// @brief This function finds the top and bottom end of a detector segment in
/// local coordinates
///
/// @param [in] local Local position of the Cluster
/// @param [in] segment Segmentation of the detector element
///
/// @return Pair containing the top and bottom end
inline std::pair<Vector2, Vector2> findLocalTopAndBottomEnd(
    const Vector2& local, const Segmentation* segment) {
  auto& binData = segment->binUtility().binningData();
  auto& boundariesX = binData[0].boundaries();
  auto& boundariesY = binData[1].boundaries();

  // Search the x-/y-bin of the Cluster
  size_t binX = binData[0].searchLocal(local);
  size_t binY = binData[1].searchLocal(local);

  // Storage of the local top (first) and bottom (second) end
  std::pair<Vector2, Vector2> topBottomLocal;

  if (boundariesX[binX + 1] - boundariesX[binX] <
      boundariesY[binY + 1] - boundariesY[binY]) {
    topBottomLocal.first = {(boundariesX[binX] + boundariesX[binX + 1]) / 2,
                            boundariesY[binY + 1]};
    topBottomLocal.second = {(boundariesX[binX] + boundariesX[binX + 1]) / 2,
                             boundariesY[binY]};
  } else {
    // Set the top and bottom end of the strip in local coordinates
    topBottomLocal.first = {boundariesX[binX],
                            (boundariesY[binY] + boundariesY[binY + 1]) / 2};
    topBottomLocal.second = {boundariesX[binX + 1],
                             (boundariesY[binY] + boundariesY[binY + 1]) / 2};
  }

  return topBottomLocal;
}

// /// @brief Calculates a space point whithout using the vertex
// /// @note This is mostly to resolve space points from cosmic data
// /// @param a vector to the top end of the first SDE
// /// @param c vector to the top end of the second SDE
// /// @param q vector from the bottom to the top end of the first SDE
// /// @param r vector from the bottom to the top end of the second SDE
// /// @return parameter that indicates the location of the space point; returns
// /// 1. if it failed
// /// @note The meaning of the parameter is explained in more detail in the
// /// function body
inline double calcPerpendicularProjection(const Vector3& a, const Vector3& c,
                                          const Vector3& q, const Vector3& r) {
  /// This approach assumes that no vertex is available. This option aims to
  /// approximate the space points from cosmic data.
  /// The underlying assumption is that the best point is given by the  closest
  /// distance between both lines describing the SDEs.
  /// The point x on the first SDE is parametrized as a + lambda0 * q with  the
  /// top end a of the strip and the vector q = a - b(ottom end of the  strip).
  /// An analogous parametrization is performed of the second SDE with y = c  +
  /// lambda1 * r.
  /// x get resolved by resolving lambda0 from the condition that |x-y| is  the
  /// shortest distance between two skew lines.

  Vector3 ac = c - a;
  double qr = q.dot(r);
  double denom = q.dot(q) - qr * qr;

  // Check for numerical stability
  if (fabs(denom) > 10e-7) {
    // Return lambda0
    return (ac.dot(r) * qr - ac.dot(q) * r.dot(r)) / denom;
  }
  // lambda0 is in the interval [-1,0]. This return serves as error check.
  return 1.;
}

/// @brief This function tests if a space point can be estimated by a more
/// tolerant treatment of construction. In fact, this function indirectly
/// allows shifts of the vertex.
///
/// @param [in] spaPoPa container that stores geometric parameters and rules of
/// the space point formation
/// @param [in] stripLengthGapTolerance Tolerance scaling factor of the gap
/// between strip detector elements
///
/// @return indicator if the test was successful
inline bool recoverSpacePoint(SpacePointParameters& spaPoPa,
                              double stripLengthGapTolerance) {
  /// Consider some cases that would allow an easy exit
  // Check if the limits are allowed to be increased
  if (stripLengthGapTolerance <= 0.) {
    return false;
  }
  spaPoPa.qmag = spaPoPa.q.norm();
  // Increase the limits. This allows a check if the point is just slightly
  // outside the SDE
  spaPoPa.limitExtended =
      spaPoPa.limit + stripLengthGapTolerance / spaPoPa.qmag;
  // Check if m is just slightly outside
  if (fabs(spaPoPa.m) > spaPoPa.limitExtended) {
    return false;
  }
  // Calculate n if not performed previously
  if (spaPoPa.n == 0.) {
    spaPoPa.n = -spaPoPa.t.dot(spaPoPa.qs) / spaPoPa.r.dot(spaPoPa.qs);
  }
  // Check if n is just slightly outside
  if (fabs(spaPoPa.n) > spaPoPa.limitExtended) {
    return false;
  }

  /// The following code considers an overshoot of m and n in the same direction
  /// of their SDE. The term "overshoot" represents the amount of m or n outside
  /// its regular interval (-1, 1).
  /// It calculates which overshoot is worse. In order to compare both, the
  /// overshoot in n is projected onto the first surface by considering the
  /// normalized projection of r onto q.
  /// This allows a rescaling of the overshoot. The worse overshoot will be set
  /// to +/-1, the parameter with less overshoot will be moved towards 0 by the
  /// worse overshoot.
  /// In order to treat both SDEs equally, the rescaling eventually needs to be
  /// performed several times. If these shifts allows m and n to be in the
  /// limits, the space point can be stored.
  /// @note This shift can be understood as a shift of the particle's
  /// trajectory. This is leads to a shift of the vertex. Since these two points
  /// are treated independently from other cluster, it is also possible to
  /// consider this as a change in the slope of the particle's trajectory.
  ///  The would also move the vertex position.

  // Calculate the scaling factor to project lengths of the second SDE on the
  // first SDE
  double secOnFirstScale =
      spaPoPa.q.dot(spaPoPa.r) / (spaPoPa.qmag * spaPoPa.qmag);
  // Check if both overshoots are in the same direction
  if (spaPoPa.m > 1. && spaPoPa.n > 1.) {
    // Calculate the overshoots
    double mOvershoot = spaPoPa.m - 1.;
    double nOvershoot =
        (spaPoPa.n - 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    spaPoPa.m -= biggerOvershoot;
    spaPoPa.n -= (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spaPoPa.m) < spaPoPa.limit && fabs(spaPoPa.n) < spaPoPa.limit;
  }
  // Check if both overshoots are in the same direction
  if (spaPoPa.m < -1. && spaPoPa.n < -1.) {
    // Calculate the overshoots
    double mOvershoot = -(spaPoPa.m + 1.);
    double nOvershoot =
        -(spaPoPa.n + 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    spaPoPa.m += biggerOvershoot;
    spaPoPa.n += (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spaPoPa.m) < spaPoPa.limit && fabs(spaPoPa.n) < spaPoPa.limit;
  }
  // No solution could be found
  return false;
}

/// @brief This function performs a straight forward calculation of a space
/// point and returns whether it was succesful or not.
///
/// @param [in] stripEnds1 Top and bottom end of the first strip detector
/// element
/// @param [in] stripEnds1 Top and bottom end of the second strip detector
/// element
/// @param [in] posVertex Position of the vertex
/// @param [in, out] spaPoPa Data container of the calculations
/// @param [in] stripLengthTolerance Tolerance scaling factor on the strip
/// detector element length
///
/// @return Boolean statement whether the space point calculation was succesful
inline bool calculateSpacePoint(const std::pair<Vector3, Vector3>& stripEnds1,
                                const std::pair<Vector3, Vector3>& stripEnds2,
                                const Vector3& posVertex,
                                SpacePointParameters& spaPoPa,
                                const double stripLengthTolerance) {
  /// The following algorithm is meant for finding the position on the first
  /// strip if there is a corresponding Cluster on the second strip. The
  /// resulting point is a point x on the first surfaces. This point is
  /// along a line between the points a (top end of the strip)
  /// and b (bottom end of the strip). The location can be parametrized as
  /// 	2 * x = (1 + m) a + (1 - m) b
  /// as function of the scalar m. m is a parameter in the interval
  /// -1 < m < 1 since the hit was on the strip. Furthermore, the vector
  /// from the vertex to the Cluster on the second strip y is needed to be a
  /// multiple k of the vector from vertex to the hit on the first strip x.
  /// As a consequence of this demand y = k * x needs to be on the
  /// connecting line between the top (c) and bottom (d) end of
  /// the second strip. If both clusters correspond to each other, the
  /// condition
  /// 	y * (c X d) = k * x (c X d) = 0 ("X" represents a cross product)
  /// needs to be fulfilled. Inserting the first equation into this
  /// equation leads to the condition for m as given in the following
  /// algorithm and therefore to the calculation of x.
  /// The same calculation can be repeated for y. Its corresponding
  /// parameter will be named n.

  spaPoPa.s = stripEnds1.first + stripEnds1.second - 2 * posVertex;
  spaPoPa.t = stripEnds2.first + stripEnds2.second - 2 * posVertex;
  spaPoPa.qs = spaPoPa.q.cross(spaPoPa.s);
  spaPoPa.rt = spaPoPa.r.cross(spaPoPa.t);
  spaPoPa.m = -spaPoPa.s.dot(spaPoPa.rt) / spaPoPa.q.dot(spaPoPa.rt);

  // Set the limit for the parameter
  if (spaPoPa.limit == 1. && stripLengthTolerance != 0.) {
    spaPoPa.limit = 1. + stripLengthTolerance;
  }

  // Check if m and n can be resolved in the interval (-1, 1)
  return (fabs(spaPoPa.m) <= spaPoPa.limit &&
          fabs(spaPoPa.n = -spaPoPa.t.dot(spaPoPa.qs) /
                           spaPoPa.r.dot(spaPoPa.qs)) <= spaPoPa.limit);
}
}  // namespace detail
}  // namespace Acts

template <typename spacepoint_t>
Acts::SpacePointBuilder<spacepoint_t>::SpacePointBuilder(
    Acts::SpacePointBuilderConfig cfg, std::unique_ptr<const Logger> logger)
    : m_config(cfg), m_logger(std::move(logger)) {}

template <typename spacepoint_t>
std::pair<Acts::Vector3, Acts::Vector2>
Acts::SpacePointBuilder<spacepoint_t>::globalCoords(
    const GeometryContext& gctx, const Measurement& meas) const {
  const auto& slink =
      std::visit([](const auto& x) { return &x.sourceLink(); }, meas);
  const auto geoId = slink->geometryId();

  const Acts::Surface* surface = m_config.trackingGeometry->findSurface(geoId);
  auto [localPos, localCov] = std::visit(
      [](const auto& measurement) {
        auto expander = measurement.expander();
        Acts::BoundVector par = expander * measurement.parameters();
        Acts::BoundSymMatrix cov =
            expander * measurement.covariance() * expander.transpose();
        // extract local position
        Acts::Vector2 lpar(par[Acts::eBoundLoc0], par[Acts::eBoundLoc1]);
        // extract local position covariance.
        Acts::SymMatrix2 lcov =
            cov.block<2, 2>(Acts::eBoundLoc0, Acts::eBoundLoc0);
        return std::make_pair(lpar, lcov);
      },
      meas);

  // transform local position to global coordinates
  Acts::Vector3 globalFakeMom(1, 1, 1);

  Acts::Vector3 globalPos =
      surface->localToGlobal(gctx, localPos, globalFakeMom);
  Acts::RotationMatrix3 rotLocalToGlobal =
      surface->referenceFrame(gctx, globalPos, globalFakeMom);

  // the space point requires only the variance of the transverse and
  // longitudinal position. reduce computations by transforming the
  // covariance directly from local to rho/z.
  //
  // compute Jacobian from global coordinates to rho/z
  //
  //         rho = sqrt(x² + y²)
  // drho/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
  //             = 2 * {x,y} / r
  //       dz/dz = 1 (duuh!)
  //
  auto x = globalPos[Acts::ePos0];
  auto y = globalPos[Acts::ePos1];
  auto scale = 2 / std::hypot(x, y);
  Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
  jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
  jacXyzToRhoZ(1, Acts::ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  Acts::ActsMatrix<2, 2> jac =
      jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
  // compute rho/z variance
  Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  auto gcov = Acts::Vector2(var[0], var[1]);
  return std::make_pair(globalPos, gcov);
}

template <typename spacepoint_t>
void Acts::SpacePointBuilder<spacepoint_t>::calculateSingleHitSpacePoints(
    const GeometryContext& gctx,
    const std::vector<const Measurement*>& measurements,
    std::vector<spacepoint_t>& spacePointStorage) const {
  for (const auto& meas : measurements) {
    auto [gPos, gCov] = globalCoords(gctx, *meas);

    std::cout << "gpos " << gPos << std::endl;

    const auto& slink =
        std::visit([](const auto& x) { return &x.sourceLink(); }, *meas);

    std::vector<const SourceLink*> slinks;
    slinks.emplace_back(slink);
    auto sp = spacepoint_t(gPos, gCov[0], gCov[1], slinks);

    spacePointStorage.emplace_back(sp);
  }
}

template <typename spacepoint_t>
void Acts::SpacePointBuilder<spacepoint_t>::calculateSpacePoints(
    const GeometryContext& gctx, std::vector<spacepoint_t>& spacePointStorage,
    const std::vector<const Measurement*>* FrontMeasurements,
    const std::vector<const Measurement*>* BackMeasurements) const {
  if (FrontMeasurements == nullptr) {
    ACTS_WARNING(" No measurements found in the SP formation");
    return;
  }

  if (BackMeasurements == nullptr) {  // pixel SP building
    calculateSingleHitSpacePoints(gctx, *FrontMeasurements, spacePointStorage);

  } else {  // strip SP building

    std::vector<std::pair<const Measurement*, const Measurement*>>
        measurementPairs;
    makeMeasurementPairs(gctx, *FrontMeasurements, *BackMeasurements,
                         measurementPairs);
    std::cout << "number of measurement pairs " << measurementPairs.size()
              << std::endl;

    calculateDoubleHitSpacePoints(gctx, measurementPairs, spacePointStorage);
  }
}

template <typename spacepoint_t>
void Acts::SpacePointBuilder<spacepoint_t>::makeMeasurementPairs(
    const GeometryContext& gctx,
    const std::vector<const Measurement*>& measurementsFront,
    const std::vector<const Measurement*>& measurementsBack,
    std::vector<std::pair<const Measurement*, const Measurement*>>&
        measurementPairs) const {
  std::cout << "make measurement pairs" << std::endl;
  // Return if no Measurements are given in a vector
  if (measurementsFront.empty() || measurementsBack.empty()) {
    return;
  }
  // Declare helper variables
  double currentDiff;
  double diffMin;
  unsigned int measurementMinDist;

  // Walk through all Measurements on both surfaces
  for (unsigned int iMeasurementsFront = 0;
       iMeasurementsFront < measurementsFront.size(); iMeasurementsFront++) {
    // Set the closest distance to the maximum of double
    diffMin = std::numeric_limits<double>::max();
    // Set the corresponding index to an element not in the list of Measurements
    measurementMinDist = measurementsBack.size();
    for (unsigned int iMeasurementsBack = 0;
         iMeasurementsBack < measurementsBack.size(); iMeasurementsBack++) {
      // auto gpos_front = globalPos(gctx,
      // measurementsFront[iMeasurementsFront]);
      auto [gposFront, gcovFront] =
          globalCoords(gctx, *(measurementsFront[iMeasurementsFront]));
      auto [gposBack, gcovBack] =
          globalCoords(gctx, *(measurementsBack[iMeasurementsBack]));

      currentDiff = detail::differenceOfMeasurementsChecked(
          gposFront, gposBack, m_config.vertex, m_config.diffDist,
          m_config.diffPhi2, m_config.diffTheta2);

      // Store the closest Measurements (distance and index) calculated so far
      if (currentDiff < diffMin && currentDiff >= 0.) {
        diffMin = currentDiff;
        measurementMinDist = iMeasurementsBack;
      }
    }

    // Store the best (=closest) result
    if (measurementMinDist < measurementsBack.size()) {
      std::pair<const Measurement*, const Measurement*> measurementPair;
      measurementPair = std::make_pair(measurementsFront[iMeasurementsFront],
                                       measurementsBack[measurementMinDist]);
      measurementPairs.push_back(measurementPair);
    }
  }
}

template <typename spacepoint_t>
void Acts::SpacePointBuilder<spacepoint_t>::calculateDoubleHitSpacePoints(
    const Acts::GeometryContext& gctx,
    const std::vector<std::pair<const Measurement*, const Measurement*>>&
        measurementPairs,
    std::vector<spacepoint_t>& spacePoints) const {
  // Source of algorithm: Athena, SiSpacePointMakerTool::makeSCT_SpacePoint()

  detail::SpacePointParameters spaPoPa;

  // Walk over every found candidate pair
  for (const auto& mp : measurementPairs) {
    // Calculate the ends of the SDEs

    const auto& ends1 = endsOfStrip(gctx, *(mp.first));

    const auto& ends2 = endsOfStrip(gctx, *(mp.second));
    std::cout << "ends obtained " << std::endl;
    spaPoPa.q = ends1.first - ends1.second;
    spaPoPa.r = ends2.first - ends2.second;

    const auto& slink_front = getSourceLink(*(mp.first));
    const auto& slink_back = getSourceLink(*(mp.second));

    std::vector<const SourceLink*> slinks;
    slinks.emplace_back(slink_front);
    slinks.emplace_back(slink_back);
    // std::vector<size_t> measurementIndices = {id_front, id_back};
    double resultPerpProj;

    if (m_config.usePerpProj) {  // for cosmic without vertex constraint
      resultPerpProj = detail::calcPerpendicularProjection(
          ends1.first, ends2.first, spaPoPa.q, spaPoPa.r);
      if (resultPerpProj <= 0.) {
        Vector3 pos = ends1.first + resultPerpProj * spaPoPa.q;
        double theta = acos(spaPoPa.q.dot(spaPoPa.r) /
                            (spaPoPa.q.norm() * spaPoPa.r.norm()));
        const auto gcov = getGlobalVars(gctx, *(mp.first), *(mp.second), theta);
        const double varRho = gcov[0];
        const double varZ = gcov[1];

        auto sp = spacepoint_t(pos, varRho, varZ, std::move(slinks));
        spacePoints.push_back(std::move(sp));
        continue;
      }
    }

    bool spFound = calculateSpacePoint(ends1, ends2, m_config.vertex, spaPoPa,
                                       m_config.stripLengthTolerance);

    if (!spFound)
      spFound =
          detail::recoverSpacePoint(spaPoPa, m_config.stripLengthGapTolerance);

    if (spFound) {
      Vector3 pos = 0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q);

      double theta = acos(spaPoPa.q.dot(spaPoPa.r) /
                          (spaPoPa.q.norm() * spaPoPa.r.norm()));

      const auto gcov = getGlobalVars(gctx, *(mp.first), *(mp.second), theta);
      const double varRho = gcov[0];
      const double varZ = gcov[1];
      auto sp = spacepoint_t(pos, varRho, varZ, std::move(slinks));
      spacePoints.push_back(std::move(sp));
    }
  }
}

template <typename spacepoint_t>
std::pair<Acts::Vector3, Acts::Vector3>
Acts::SpacePointBuilder<spacepoint_t>::endsOfStrip(
    const Acts::GeometryContext& gctx, const Measurement& measurement) const {
  // Calculate the local coordinates of the Cluster

  Acts::Vector3 topGlobal(0, 0, 0);
  Acts::Vector3 bottomGlobal(0, 0, 0);

  auto localPos = getLocalPos(measurement);

  const auto& slink = getSourceLink(measurement);
  const auto geoId = slink->geometryId();
  const Acts::Surface* surface = m_config.trackingGeometry->findSurface(geoId);

  auto detectorElement = dynamic_cast<const Acts::IdentifiedDetectorElement*>(
      surface->associatedDetectorElement());

  if (!detectorElement && detectorElement->digitizationModule()) {
    ACTS_ERROR(" No detector element found for the strip SP formation");

    return std::make_pair(topGlobal, bottomGlobal);
  }

  auto digitationModule = detectorElement->digitizationModule();
  const Acts::Segmentation& segment = digitationModule->segmentation();

  std::pair<Acts::Vector2, Acts::Vector2> topBottomLocal =
      detail::findLocalTopAndBottomEnd(localPos, &segment);

  // Calculate the global coordinates of the top and bottom end of the strip
  Acts::Vector3 globalFakeMom(1, 1, 1);
  topGlobal = surface->localToGlobal(gctx, topBottomLocal.first, globalFakeMom);
  bottomGlobal =
      surface->localToGlobal(gctx, topBottomLocal.second, globalFakeMom);

  // Return the top and bottom end of the strip in global coordinates
  return std::make_pair(topGlobal, bottomGlobal);
}

template <typename spacepoint_t>
Acts::Vector2 Acts::SpacePointBuilder<spacepoint_t>::getGlobalVars(
    const Acts::GeometryContext& gctx, const Measurement& meas_front,
    const Measurement& meas_back, const double theta) const {
  std::cout << "getGlobalVars" << std::endl;
  const auto var1 = getLocVar(meas_front);
  const auto var2 = getLocVar(meas_back);
  // strip1 and strip2 are tilted at +/- theta/2

  double sigma_x = std::hypot(var1, var2) / (2 * sin(theta * 0.5));
  double sigma_y = std::hypot(var1, var2) / (2 * cos(theta * 0.5));

  // projection to the surface with strip1.
  double sig_x1 = sigma_x * cos(0.5 * theta) + sigma_y * sin(0.5 * theta);
  double sig_y1 = sigma_y * cos(0.5 * theta) + sigma_x * sin(0.5 * theta);
  Acts::SymMatrix2 lcov;
  lcov << sig_x1, 0, 0, sig_y1;

  auto [localPos, localCov] = getLocalPosCov(meas_front);

  const auto& slink_meas1 =
      std::visit([](const auto& x) { return &x.sourceLink(); }, meas_front);

  const auto geoId = slink_meas1->geometryId();

  auto gcov = globalCov(gctx, geoId, localPos, lcov);

  return gcov;
}

template <typename spacepoint_t>
double Acts::SpacePointBuilder<spacepoint_t>::getLocVar(
    const Measurement& meas) const {
  std::cout << "getLocalVar" << std::endl;
  auto cov = std::visit(
      [](const auto& x) {
        auto expander = x.expander();
        Acts::BoundSymMatrix bcov =
            expander * x.covariance() * expander.transpose();
        Acts::SymMatrix2 lcov =
            bcov.block<2, 2>(Acts::eBoundLoc0, Acts::eBoundLoc0);
        return lcov;
      },
      meas);
  return cov(0, 0);
}

template <typename spacepoint_t>
Acts::Vector2 Acts::SpacePointBuilder<spacepoint_t>::getLocalPos(
    const Measurement& meas) const {
  std::cout << "getLocalPos" << std::endl;
  Acts::Vector2 localPos = std::visit(
      [](const auto& x) {
        auto expander = x.expander();
        Acts::BoundVector par = expander * x.parameters();
        std::cout << par << std::endl;
        Acts::Vector2 local(par[Acts::BoundIndices::eBoundLoc0],
                            par[Acts::BoundIndices::eBoundLoc1]);
        return local;
      },
      meas);
  return localPos;
}

template <typename spacepoint_t>
std::pair<Acts::Vector2, Acts::SymMatrix2>
Acts::SpacePointBuilder<spacepoint_t>::getLocalPosCov(
    const Measurement& meas) const {
  std::cout << "getLocalPosCov" << std::endl;
  return std::visit(
      [](const auto& x) {
        auto expander = x.expander();
        Acts::BoundVector par = expander * x.parameters();
        Acts::BoundSymMatrix cov =
            expander * x.covariance() * expander.transpose();
        Acts::Vector2 lpar(par[Acts::BoundIndices::eBoundLoc0],
                           par[Acts::BoundIndices::eBoundLoc1]);
        Acts::SymMatrix2 lcov =
            cov.block<2, 2>(Acts::eBoundLoc0, Acts::eBoundLoc0);
        return std::make_pair(lpar, lcov);
      },
      meas);
}

template <typename spacepoint_t>
Acts::Vector2 Acts::SpacePointBuilder<spacepoint_t>::globalCov(
    const Acts::GeometryContext& gctx, const Acts::GeometryIdentifier& geoId,
    const Acts::Vector2& localPos, const Acts::SymMatrix2& localCov) const {
  std::cout << "globalCov" << std::endl;

  Acts::Vector3 globalFakeMom(1, 1, 1);

  const Acts::Surface* surface = m_config.trackingGeometry->findSurface(geoId);

  Acts::Vector3 globalPos =
      surface->localToGlobal(gctx, localPos, globalFakeMom);
  Acts::RotationMatrix3 rotLocalToGlobal =
      surface->referenceFrame(gctx, globalPos, globalFakeMom);

  auto x = globalPos[Acts::ePos0];
  auto y = globalPos[Acts::ePos1];
  auto scale = 2 / std::hypot(x, y);
  Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
  jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
  jacXyzToRhoZ(1, Acts::ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  Acts::ActsMatrix<2, 2> jac =
      jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
  // compute rho/z variance
  Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  auto gcov = Acts::Vector2(var[0], var[1]);

  return gcov;
}

template <typename spacepoint_t>
const Acts::SourceLink* Acts::SpacePointBuilder<spacepoint_t>::getSourceLink(
    const Measurement meas) const {
  const Acts::SourceLink* slink =
      std::visit([](const auto& x) { return &x.sourceLink(); }, meas);
  return slink;
}
