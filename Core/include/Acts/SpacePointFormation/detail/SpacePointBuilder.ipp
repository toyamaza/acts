// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
namespace Acts {


Result<void> testRecoverSpacePoint(
    SpacePointParameters& spParams, double stripLengthGapTolerance)  {
  std::error_code m_error;
  
  // Consider some cases that would allow an easy exit
  // Check if the limits are allowed to be increased
  std::cout << "recovering space point " << std::endl;

  if (stripLengthGapTolerance <= 0.) {
    std::cout << "failure 1 (tolerance <=0 )" << std::endl;
    return Result<void>::failure(m_error);
  }

  spParams.mag_firstBtmToTop = spParams.firstBtmToTop.norm();
  // Increase the limits. This allows a check if the point is just slightly
  // outside the SDE
  spParams.limitExtended =
      spParams.limit + stripLengthGapTolerance / spParams.mag_firstBtmToTop;
  std::cout << "limExtended " << spParams.limitExtended << std::endl;
  // Check if m is just slightly outside
  if (fabs(spParams.m) > spParams.limitExtended) {
    std::cout << "failure 2" << std::endl;    
    return Result<void>::failure(m_error);
  }
  // Calculate n if not performed previously
  if (spParams.n == 0.) {
    spParams.n =
        -spParams.vtxToSecondMid2.dot(spParams.firstBtmToTopXvtxToFirstMid2) /
        spParams.secondBtmToTop.dot(spParams.firstBtmToTopXvtxToFirstMid2);
  }

  // Check if n is just slightly outside
  if (fabs(spParams.n) > spParams.limitExtended) {
    std::cout << "failure 3 (n < limExt )" << std::endl;    
    return Result<void>::failure(m_error);
  }
  
  double secOnFirstScale =
      spParams.firstBtmToTop.dot(spParams.secondBtmToTop) /
      (spParams.mag_firstBtmToTop * spParams.mag_firstBtmToTop);
  // Check if both overshoots are in the same direction
  if (spParams.m > 1. && spParams.n > 1.) {
    // Calculate the overshoots
    double mOvershoot = spParams.m - 1.;
    double nOvershoot =
        (spParams.n - 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    spParams.m -= biggerOvershoot;
    spParams.n -= (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point

    if (fabs(spParams.m) < spParams.limit &&
        fabs(spParams.n) < spParams.limit) {
      return Result<void>::success();
    } else {
    std::cout << "failure 4" << std::endl;      
      return Result<void>::failure(m_error);
    }
  }
  // Check if both overshoots are in the same direction
  if (spParams.m < -1. && spParams.n < -1.) {
    // Calculate the overshoots
    double mOvershoot = -(spParams.m + 1.);
    double nOvershoot =
        -(spParams.n + 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    spParams.m += biggerOvershoot;
    spParams.n += (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    if (fabs(spParams.m) < spParams.limit &&
        fabs(spParams.n) < spParams.limit) {
      return Result<void>::success();
    }
  }
  // No solution could be found
    std::cout << "failure 5" << std::endl;  
  return Result<void>::failure(m_error);
}


  
template <typename spacepoint_t>
SpacePointBuilder<spacepoint_t>::SpacePointBuilder(
    SpacePointBuilderConfig cfg,
    std::function<
        spacepoint_t(Acts::Vector3, Acts::Vector2,
                     boost::container::static_vector<const SourceLink*, 2>)>
        func,
    std::unique_ptr<const Logger> logger)
    : m_config(cfg), m_spConstructor(func), m_logger(std::move(logger)) {
  m_spUtility = std::make_shared<SpacePointUtility>(cfg);
}

template <typename spacepoint_t>
template <template <typename...> typename container_t>
void SpacePointBuilder<spacepoint_t>::buildSpacePoint(
    const GeometryContext& gctx,
    const std::vector<const Measurement*>& measurements,
    const SpacePointBuilderOptions& opt,
    std::back_insert_iterator<container_t<spacepoint_t>> spacePointIt) const {
  const unsigned int num_meas = measurements.size();

  Acts::Vector3 gPos = Acts::Vector3::Zero();
  Acts::Vector2 gCov = Acts::Vector2::Zero();

  if (num_meas == 1) {  // pixel SP formation

    auto gPosCov = m_spUtility->globalCoords(gctx, *(measurements[0]));
    gPos = gPosCov.first;
    gCov = gPosCov.second;

  } else if (num_meas == 2) {  // strip SP formation

    const auto& ends1 = opt.stripEndsPair.first;
    const auto& ends2 = opt.stripEndsPair.second;

    Acts::SpacePointParameters spParams;

    if (!m_config.usePerpProj) {  // default strip SP building

      auto spFound = m_spUtility->calculateStripSPPosition(
          ends1, ends2, opt.vertex, spParams, opt.stripLengthTolerance);

      if (!spFound.ok()) {
        ACTS_VERBOSE(
            "SP formation: First attempt failed. Trying to recover SP");
        std::cout << "m, n, limit " << fabs(spParams.m) << " "
                  << fabs(spParams.n) << " " << spParams.limit << std::endl;

        spFound = m_spUtility->recoverSpacePoint(spParams,
                                                 opt.stripLengthGapTolerance);
        testRecoverSpacePoint(spParams,
			      opt.stripLengthGapTolerance);

        // ------------------------ tmp

        // ------------------------
      }

      if (!spFound.ok()) {
        ACTS_VERBOSE("SP formation: no SP created for this pair")
        return;
      }

      gPos = 0.5 *
             (ends1.first + ends1.second + spParams.m * spParams.firstBtmToTop);

    } else {  // for cosmic without vertex constraint

      auto resultPerpProj =
          m_spUtility->calcPerpendicularProjection(ends1, ends2, spParams);

      if (!resultPerpProj.ok())
        return;
      gPos = ends1.first + resultPerpProj.value() * spParams.firstBtmToTop;
    }

    double theta =
        acos(spParams.firstBtmToTop.dot(spParams.secondBtmToTop) /
             (spParams.firstBtmToTop.norm() * spParams.secondBtmToTop.norm()));

    gCov = m_spUtility->calcRhoZVars(gctx, *(measurements.at(0)),
                                     *(measurements.at(1)), gPos, theta);

  } else {
    ACTS_ERROR("More than 2 measurements are given for a space point.");
  }

  boost::container::static_vector<const SourceLink*, 2> slinks;
  for (const auto& meas : measurements) {
    const auto& slink =
        std::visit([](const auto& x) { return &x.sourceLink(); }, *meas);
    slinks.emplace_back(slink);
  }

  spacePointIt = m_spConstructor(gPos, gCov, std::move(slinks));
}

template <typename spacepoint_t>
void SpacePointBuilder<spacepoint_t>::makeMeasurementPairs(
    const GeometryContext& gctx, const StripPairOptions& pairOpt,
    const std::vector<const Measurement*>& measurementsFront,
    const std::vector<const Measurement*>& measurementsBack,
    std::vector<std::pair<const Measurement*, const Measurement*>>&
        measurementPairs) const {
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
      auto [gposFront, gcovFront] = m_spUtility->globalCoords(
          gctx, *(measurementsFront[iMeasurementsFront]));
      auto [gposBack, gcovBack] = m_spUtility->globalCoords(
          gctx, *(measurementsBack[iMeasurementsBack]));

      auto res = m_spUtility->differenceOfMeasurementsChecked(
          gposFront, gposBack, pairOpt.vertex, pairOpt.diffDist,
          pairOpt.diffPhi2, pairOpt.diffTheta2);
      if (!res.ok())
        continue;

      currentDiff = res.value();

      // Store the closest Measurements (distance and index) calculated so far
      if (currentDiff < diffMin && currentDiff >= 0.) {
        diffMin = currentDiff;
        measurementMinDist = iMeasurementsBack;
      }
    }

    // Store the best (=closest) result
    if (measurementMinDist < measurementsBack.size()) {
      std::pair<const Measurement*, const Measurement*> measurementPair =
          std::make_pair(measurementsFront[iMeasurementsFront],
                         measurementsBack[measurementMinDist]);
      measurementPairs.push_back(measurementPair);
    }
  }
}

}  // namespace Acts
