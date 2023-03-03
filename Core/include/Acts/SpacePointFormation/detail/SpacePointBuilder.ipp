// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
namespace Acts {

template <typename spacepoint_t>
Result<void> SpacePointBuilder<spacepoint_t>::testRecoverSpacePoint(
    SpacePointParameters& spParams, double stripLengthGapTolerance) const {
  std::error_code m_error;

  // Consider some cases that would allow an easy exit
  // Check if the limits are allowed to be increased

  if (stripLengthGapTolerance <= 0.) {
    return Result<void>::failure(m_error);
  }

  spParams.mag_firstBtmToTop = spParams.firstBtmToTop.norm();
  // Increase the limits. This allows a check if the point is just slightly
  // outside the SDE
  spParams.limitExtended =
      spParams.limit + stripLengthGapTolerance / spParams.mag_firstBtmToTop;
  std::cout << "limExtended =  (lim + tolerance)/ mag  = "
            << spParams.limitExtended << " =  " << spParams.limit << " + "
            << stripLengthGapTolerance << " / " << spParams.mag_firstBtmToTop
            << std::endl;
  // Check if m is just slightly outside
  if (fabs(spParams.m) > spParams.limitExtended) {
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
  return Result<void>::failure(m_error);
}

template <typename spacepoint_t>
SpacePointBuilder<spacepoint_t>::SpacePointBuilder(
    const SpacePointBuilderConfig& cfg,
    std::function<spacepoint_t(Acts::Vector3, Acts::Vector2,
                               boost::container::static_vector<SourceLink, 2>)>
        func,
    std::unique_ptr<const Logger> logger)
    : m_config(cfg), m_spConstructor(func), m_logger(std::move(logger)) {
  m_spUtility = std::make_shared<SpacePointUtility>(cfg);
}

template <typename spacepoint_t>
template <template <typename...> typename container_t>
void SpacePointBuilder<spacepoint_t>::buildSpacePoint(
    const GeometryContext& gctx, const std::vector<SourceLink>& sourceLinks,
    const SpacePointBuilderOptions& opt,
    std::back_insert_iterator<container_t<spacepoint_t>> spacePointIt) const {
  const unsigned int num_slinks = sourceLinks.size();

  Acts::Vector3 gPos = Acts::Vector3::Zero();
  Acts::Vector2 gCov = Acts::Vector2::Zero();

  if (num_slinks == 1) {  // pixel SP formation
    auto slink = sourceLinks.at(0);
    auto [param, cov] = opt.paramCovAccessor(sourceLinks.at(0));
    auto gPosCov = m_spUtility->globalCoords(gctx, slink, param, cov);
    gPos = gPosCov.first;
    gCov = gPosCov.second;
  } else if (num_slinks == 2) {  // strip SP formation

    const auto& ends1 = opt.stripEndsPair.first;
    const auto& ends2 = opt.stripEndsPair.second;

    Acts::SpacePointParameters spParams;

    if (!m_config.usePerpProj) {  // default strip SP building

      auto spFound = m_spUtility->calculateStripSPPosition(
          ends1, ends2, opt.vertex, spParams, opt.stripLengthTolerance);

      if (!spFound.ok()) {
        ACTS_VERBOSE(
            "SP formation: First attempt failed. Trying to recover SP");

        spFound = m_spUtility->recoverSpacePoint(spParams,
                                                 opt.stripLengthGapTolerance);
        testRecoverSpacePoint(spParams, opt.stripLengthGapTolerance);

        // ------------------------ tmp

        // ------------------------
      }

      if (!spFound.ok()) {
        return;
      }

      gPos = 0.5 *
             (ends1.first + ends1.second + spParams.m * spParams.firstBtmToTop);

    } else {  // for cosmic without vertex constraint

      auto resultPerpProj =
          m_spUtility->calcPerpendicularProjection(ends1, ends2, spParams);

      if (!resultPerpProj.ok()) {
        return;
      }
      gPos = ends1.first + resultPerpProj.value() * spParams.firstBtmToTop;
    }

    double theta =
        acos(spParams.firstBtmToTop.dot(spParams.secondBtmToTop) /
             (spParams.firstBtmToTop.norm() * spParams.secondBtmToTop.norm()));

    gCov = m_spUtility->calcRhoZVars(gctx, sourceLinks.at(0), sourceLinks.at(1),
                                     opt.paramCovAccessor, gPos, theta);

  } else {
    ACTS_ERROR("More than 2 sourceLinks are given for a space point.");
  }
  boost::container::static_vector<SourceLink, 2> slinks(sourceLinks.begin(),
                                                        sourceLinks.end());

  spacePointIt = m_spConstructor(gPos, gCov, std::move(slinks));
}

template <typename spacepoint_t>
void SpacePointBuilder<spacepoint_t>::makeSlinkPairs(
    const GeometryContext& gctx, const std::vector<SourceLink>& slinksFront,
    const std::vector<SourceLink>& slinksBack,
    std::vector<std::pair<SourceLink, SourceLink>>& slinkPairs,
    const StripPairOptions& pairOpt) const {
  if (slinksFront.empty() || slinksBack.empty()) {
    return;
  }
  double minDistance = 0;
  unsigned int closestIndex = 0;

  for (unsigned int i = 0; i < slinksFront.size(); i++) {
    const auto& slinkFront = slinksFront[i];
    minDistance = std::numeric_limits<double>::max();
    closestIndex = slinksBack.size();
    for (unsigned int j = 0; j < slinksBack.size(); j++) {
      const auto& slinkBack = slinksBack[j];

      const auto [paramFront, covFront] = pairOpt.paramCovAccessor(slinkFront);
      const auto [gposFront, gcovFront] =
          m_spUtility->globalCoords(gctx, slinkFront, paramFront, covFront);

      const auto [paramBack, covBack] = pairOpt.paramCovAccessor(slinkBack);
      const auto [gposBack, gcovBack] =
          m_spUtility->globalCoords(gctx, slinkBack, paramBack, covBack);

      auto res = m_spUtility->differenceOfMeasurementsChecked(
          gposFront, gposBack, pairOpt.vertex, pairOpt.diffDist,
          pairOpt.diffPhi2, pairOpt.diffTheta2);
      if (!res.ok())
        continue;
      const auto distance = res.value();
      if (distance >= 0. && distance < minDistance) {
        minDistance = distance;
        closestIndex = j;
      }
    }
    if (closestIndex < slinksBack.size()) {
      slinkPairs.emplace_back(slinksFront[i], slinksBack[closestIndex]);
    }
  }
}

}  // namespace Acts
