// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {
  
  template< typename T>
  [[deprecated]] void typecheck(T&&){}

template <typename spacepoint_t>
SpacePointBuilder<spacepoint_t>::SpacePointBuilder(
    SpacePointBuilderConfig cfg,
    std::function<
        spacepoint_t(Acts::Vector3, Acts::Vector2,
                     boost::container::static_vector<const SourceLink*, 2>)>
        func,
    std::shared_ptr<const Logger> logger)
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
      ACTS_VERBOSE("One measurement found. Perform pixel SP formation");
      std::cout << "one measurement. Pixel sp building" << std::endl;
      auto gPosCov = m_spUtility->globalCoords(gctx, *(measurements[0]));
      gPos = gPosCov.first;
      gCov = gPosCov.second;

    } else if (num_meas == 2) {  // strip SP formation
      ACTS_VERBOSE("Two measurement found. Perform strip SP formation");
      std::cout << "two measurement. Strip SP building" << std::endl;    
      const auto& ends1 = opt.stripEndsPair.first;
      const auto& ends2 = opt.stripEndsPair.second;

      // std::cout << "ends1.1 " << std::endl << ends1.first << std::endl;
      // std::cout << "ends1.2 " << std::endl << ends1.second << std::endl;    
      // std::cout << "ends2.1 " << std::endl << ends2.first << std::endl;
      // std::cout << "ends2.2 " << std::endl << ends2.second << std::endl;    


      Acts::SpacePointParameters spParams;

      if (!m_config.usePerpProj) {  // default strip SP building
	std::cout << "default strip sp building" << std::endl;

	auto spFound = m_spUtility->calculateStripSPPosition(
							     ends1, ends2, m_config.vertex, spParams,
							     m_config.stripLengthTolerance);

	std::cout << "spFound check 1 :" << spFound.ok() << std::endl;

	if (!spFound.ok()) {
	  spFound = m_spUtility->recoverSpacePoint(
						   spParams, m_config.stripLengthGapTolerance);
	}
	std::cout << "spFound check 2 :" << spFound.ok() << std::endl;

	if (!spFound.ok())
	  return;

	gPos = 0.5 *
	  (ends1.first + ends1.second + spParams.m * spParams.firstBtmToTop);

	std::cout << "gPos " << std::endl << gPos << std::endl;

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
      // std::cout << "theta " << theta << std::endl;
      // std::cout << "check0 " << std::endl;    

      auto meas0 = measurements.at(0);

      const Measurement& meas00 =  *(measurements.at(0));
      // std::cout << "check1 " << std::endl;    
      // std::cout << "check2 " << std::endl;        
      auto meas1 = measurements.at(1);
      const auto& slink_meas1 =
	std::visit([](const auto& x) { return &x.sourceLink(); }, meas00);
      std::cout << "getting geoID " << std::endl;
      const auto geoId = slink_meas1->geometryId();
      // std::cout << "getting surface " << std::endl;
      
      // std::cout << "testcov " << std::endl;
      // std::cout << std::hex << &m_config << std::endl;
      // std::cout << typeid(*(measurements[0])).name() << std::endl;

      // std::cout << "is null " <<  (m_spUtility.get() == nullptr) << std::endl;
      //      auto ttt  = m_spUtility->recoverSpacePoint( spParams, m_config.stripLengthGapTolerance);
      // std::cout << "check2.1" << std::endl;

      // const GeometryContext& geoctx = gctx;
      // std::cout << "check3" << std::endl;
      // typecheck(gctx);
      // typecheck(m_spUtility);
      // typecheck(*(measurements[0]));
      // typecheck( m_spUtility->globalCoords_tmp(gctx, *(measurements[0])) );
      std::cout << "check4.11" << std::endl;
      const Surface* surface = m_config.trackingGeometry->findSurface(geoId);
      std::cout << "check4.12" << std::endl;      
      Vector3 globalFakeMom(1, 1, 1);
      std::cout << "check4.13" << std::endl;      
      RotationMatrix3 rotLocalToGlobal =
	surface->referenceFrame(gctx, gPos, globalFakeMom);
      std::cout << "check4.14" << std::endl;      
      
      // m_spUtility->test_empty();
      // std::cout << "check4" << std::endl;
      // m_spUtility->test_meas(*(measurements[0]));
      // std::cout << "check5" << std::endl;
      // m_spUtility->test_gctx(gctx);
      // std::cout << "check6" << std::endl;      
      // m_spUtility->test_gctx_meas(gctx, *(measurements[0]));
      std::cout << "check7" << std::endl;            
    //   // auto tmpCov = m_spUtility->globalCoords_tmp(gctx, *(measurements[0]));
    //   const Measurement& measFront = *(measurements.at(0));
    //   const Measurement& measBack = *(measurements.at(1));
    //   const Vector3& globalPos = gPos;
    // std::cout << "getloc0var1" << std::endl;
    // const auto var1 = getLoc0Var(measFront);
    // std::cout << "getloc0var2" << std::endl;  
    // const auto var2 = getLoc0Var(measBack);
    // // strip1 and strip2 are tilted at +/- theta/2

    // double sigma_x = std::hypot(var1, var2) / (2 * sin(theta * 0.5));
    // double sigma_y = std::hypot(var1, var2) / (2 * cos(theta * 0.5));

    // // projection to the surface with strip1.
    // double sig_x1 = sigma_x * cos(0.5 * theta) + sigma_y * sin(0.5 * theta);
    // double sig_y1 = sigma_y * cos(0.5 * theta) + sigma_x * sin(0.5 * theta);
    // SymMatrix2 lcov;
    // lcov << sig_x1, 0, 0, sig_y1;
    // std::cout << "getting slink " << std::endl;
    // const auto& slink_meas1 =
    //   std::visit([](const auto& x) { return &x.sourceLink(); }, measFront);
    // std::cout << "getting geoID " << std::endl;
    // const auto geoId = slink_meas1->geometryId();
    // std::cout << "getting gcov " << std::endl;
    // auto gcov = rhoZCovariance(gctx, geoId, globalPos, lcov);





      
      std::cout << "calcRhoZVars" << std::endl;
      gCov = m_spUtility->calcRhoZVars(gctx, *(measurements.at(0)),
				       *(measurements.at(1)), gPos, theta);
      std::cout << "gCov " << std::endl << gCov << std::endl;    

    } else {
      ACTS_ERROR("More than 2 measurements are given for a space point.");
    }

    boost::container::static_vector<const SourceLink*, 2> slinks;
    for (const auto& meas : measurements) {
      const auto& slink =
        std::visit([](const auto& x) { return &x.sourceLink(); }, *meas);
      slinks.emplace_back(slink);
    }
    std::cout << "slink added " << std::endl;
  
    spacePointIt = m_spConstructor(gPos, gCov, std::move(slinks));
  }

template <typename spacepoint_t>
void SpacePointBuilder<spacepoint_t>::makeMeasurementPairs(
    const GeometryContext& gctx,
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
          gposFront, gposBack, m_config.vertex, m_config.diffDist,
          m_config.diffPhi2, m_config.diffTheta2);
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
