// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/TrackFitting/detail/KalmanGlobalCovariance.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

using namespace Acts::UnitLiterals;
using namespace Acts::Test;

/// Find outliers using plain distance for testing purposes.
///
/// In a real setup, the outlier classification can be much more involved, e.g.
/// by computing the weighted distance/ local chi2 value. Here, the purpose is
/// to test that the basic principle works using simplified, synthetic data.
/// Thus, the simplest possible implementation should do.
struct TestOutlierFinder {
  double distanceMax = std::numeric_limits<double>::max();

  /// Classify a measurement as a valid one or an outlier.
  ///
  /// @tparam track_state_t Type of the track state
  /// @param state The track state to classify
  /// @retval False if the measurement is not an outlier
  /// @retval True if the measurement is an outlier
  bool operator()(Acts::MultiTrajectory::ConstTrackStateProxy state) const {
    // can't determine an outlier w/o a measurement or predicted parameters
    if (not state.hasCalibrated() or not state.hasPredicted()) {
      return false;
    }
    auto residuals = state.calibrated() - state.projector() * state.predicted();
    auto distance = residuals.norm();
    return (distanceMax <= distance);
  }
};

/// Determine if the smoothing of a track should be done with or without reverse
/// filtering
struct TestReverseFilteringLogic {
  double momentumMax = std::numeric_limits<double>::max();

  /// Classify a measurement as a valid one or an outlier.
  ///
  /// @param trackState The trackState of the last measurement
  /// @retval False if we don't use the reverse filtering for the smoothing of the track
  /// @retval True if we use the reverse filtering for the smoothing of the track
  bool operator()(Acts::MultiTrajectory::ConstTrackStateProxy state) const {
    // can't determine an outlier w/o a measurement or predicted parameters
    auto momentum = fabs(1 / state.filtered()[Acts::eBoundQOverP]);
    std::cout << "momentum : " << momentum << std::endl;
    return (momentum <= momentumMax);
  }
};

// Construct a straight-line propagator.
auto makeStraightPropagator(std::shared_ptr<const Acts::TrackingGeometry> geo) {
  Acts::Navigator::Config cfg{geo};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Acts::StraightLineStepper stepper;
  return Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>(
      std::move(stepper), std::move(navigator));
}

// Construct a propagator using a constant magnetic field along z.
template <typename stepper_t>
auto makeConstantFieldPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo, double bz) {
  Acts::Navigator::Config cfg{geo};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  auto field =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, bz));
  stepper_t stepper(std::move(field));
  return Acts::Propagator<decltype(stepper), Acts::Navigator>(
      std::move(stepper), std::move(navigator));
}

// Put all this in a struct to avoid that all these objects are exposed as
// global objects in the header
struct FitterTester {
  using Rng = std::default_random_engine;

  // Context objects
  Acts::GeometryContext geoCtx;
  Acts::MagneticFieldContext magCtx;
  Acts::CalibrationContext calCtx;

  // detector geometry
  CubicTrackingGeometry geometryStore{geoCtx};
  std::shared_ptr<const Acts::TrackingGeometry> geometry = geometryStore();

  // expected number of measurements for the given detector
  constexpr static size_t nMeasurements = 6u;

  // detector resolutions
  MeasurementResolution resPixel = {MeasurementType::eLoc01, {25_um, 50_um}};
  MeasurementResolution resStrip0 = {MeasurementType::eLoc0, {100_um}};
  MeasurementResolution resStrip1 = {MeasurementType::eLoc1, {150_um}};
  MeasurementResolutionMap resolutions = {
      {Acts::GeometryIdentifier().setVolume(2), resPixel},
      {Acts::GeometryIdentifier().setVolume(3).setLayer(2), resStrip0},
      {Acts::GeometryIdentifier().setVolume(3).setLayer(4), resStrip1},
      {Acts::GeometryIdentifier().setVolume(3).setLayer(6), resStrip0},
      {Acts::GeometryIdentifier().setVolume(3).setLayer(8), resStrip1},
  };

  // simulation propagator
  Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator> simPropagator =
      makeStraightPropagator(geometry);

  //////////////////////////
  // The testing functions
  //////////////////////////

  template <typename fitter_t, typename fitter_options_t, typename parameters_t>
  void test_ZeroFieldNoSurfaceForward(const fitter_t& fitter,
                                      fitter_options_t options,
                                      const parameters_t& start, Rng& rng,
                                      const bool expected_reversed,
                                      const bool expected_smoothed) const {
    auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                           resolutions, rng);
    const auto& sourceLinks = measurements.sourceLinks;
    BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

    // this is the default option. set anyways for consistency
    options.referenceSurface = nullptr;

    auto res =
        fitter.fit(sourceLinks.begin(), sourceLinks.end(), start, options);
    BOOST_REQUIRE(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.lastMeasurementIndex, SIZE_MAX);
    BOOST_CHECK(not val.fittedParameters);
    BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
    // check the output status flags
    BOOST_CHECK(val.smoothed == expected_smoothed);
    BOOST_CHECK(val.reversed == expected_reversed);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
  }

  template <typename fitter_t, typename fitter_options_t, typename parameters_t>
  void test_ZeroFieldWithSurfaceForward(const fitter_t& fitter,
                                        fitter_options_t options,
                                        const parameters_t& start, Rng& rng,
                                        const bool expected_reversed,
                                        const bool expected_smoothed) const {
    auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                           resolutions, rng);
    const auto& sourceLinks = measurements.sourceLinks;
    BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

    // initial fitter options configured for backward filtereing mode
    // backward filtering requires a reference surface
    options.referenceSurface = &start.referenceSurface();
    // this is the default option. set anyways for consistency
    options.propagatorPlainOptions.direction = Acts::forward;

    auto res =
        fitter.fit(sourceLinks.begin(), sourceLinks.end(), start, options);
    BOOST_REQUIRE(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.lastMeasurementIndex, SIZE_MAX);
    BOOST_CHECK(val.fittedParameters);
    // check the output status flags
    BOOST_CHECK(val.smoothed == expected_smoothed);
    BOOST_CHECK(val.reversed == expected_reversed);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
    // count the number of `smoothed` states
    if (expected_reversed) {
      size_t nSmoothed = 0;
      val.fittedStates.visitBackwards(val.lastMeasurementIndex,
                                      [&nSmoothed](const auto& state) {
                                        nSmoothed += state.hasSmoothed();
                                      });
      BOOST_CHECK_EQUAL(nSmoothed, sourceLinks.size());
    }
  }

  template <typename fitter_t, typename fitter_options_t, typename parameters_t>
  void test_ZeroFieldWithSurfaceBackward(const fitter_t& fitter,
                                         fitter_options_t options,
                                         const parameters_t& start, Rng& rng,
                                         const bool expected_reversed,
                                         const bool expected_smoothed) const {
    auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                           resolutions, rng);
    const auto& sourceLinks = measurements.sourceLinks;
    BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

    // create a track near the tracker exit for outward->inward filtering
    Acts::Vector4 posOuter = start.fourPosition(geoCtx);
    posOuter[Acts::ePos0] = 3_m;
    Acts::CurvilinearTrackParameters startOuter(
        posOuter, start.unitDirection(), start.absoluteMomentum(),
        start.charge(), start.covariance());

    options.referenceSurface = &startOuter.referenceSurface();
    options.propagatorPlainOptions.direction = Acts::backward;

    auto res =
        fitter.fit(sourceLinks.begin(), sourceLinks.end(), startOuter, options);
    BOOST_CHECK(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.lastMeasurementIndex, SIZE_MAX);
    BOOST_CHECK(val.fittedParameters);
    BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
    // check the output status flags
    BOOST_CHECK(val.smoothed == expected_smoothed);
    BOOST_CHECK(val.reversed == expected_reversed);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
    // count the number of `smoothed` states
    if (expected_reversed) {
      size_t nSmoothed = 0;
      val.fittedStates.visitBackwards(val.lastMeasurementIndex,
                                      [&nSmoothed](const auto& state) {
                                        nSmoothed += state.hasSmoothed();
                                      });
      BOOST_CHECK_EQUAL(nSmoothed, sourceLinks.size());
    }
  }

  template <typename fitter_t, typename fitter_options_t, typename parameters_t>
  void test_ZeroFieldWithSurfaceAtExit(const fitter_t& fitter,
                                       fitter_options_t options,
                                       const parameters_t& start, Rng& rng,
                                       const bool expected_reversed,
                                       const bool expected_smoothed) const {
    auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                           resolutions, rng);
    const auto& sourceLinks = measurements.sourceLinks;
    BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

    // create a boundless target surface near the tracker exit
    Acts::Vector3 center(3._m, 0., 0.);
    Acts::Vector3 normal(1., 0., 0.);
    auto targetSurface =
        Acts::Surface::makeShared<Acts::PlaneSurface>(center, normal);

    options.referenceSurface = targetSurface.get();

    auto res =
        fitter.fit(sourceLinks.begin(), sourceLinks.end(), start, options);
    BOOST_REQUIRE(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.lastMeasurementIndex, SIZE_MAX);
    BOOST_CHECK(val.fittedParameters);
    BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
    // check the output status flags
    BOOST_CHECK(val.smoothed == expected_smoothed);
    BOOST_CHECK(val.reversed == expected_reversed);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
  }

  template <typename fitter_t, typename fitter_options_t, typename parameters_t>
  void test_ZeroFieldShuffled(const fitter_t& fitter, fitter_options_t options,
                              const parameters_t& start, Rng& rng,
                              const bool expected_reversed,
                              const bool expected_smoothed) const {
    auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                           resolutions, rng);
    const auto& sourceLinks = measurements.sourceLinks;
    BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

    options.referenceSurface = &start.referenceSurface();

    Acts::BoundVector parameters = Acts::BoundVector::Zero();

    // fit w/ all hits in order
    {
      auto res =
          fitter.fit(sourceLinks.begin(), sourceLinks.end(), start, options);
      BOOST_REQUIRE(res.ok());

      const auto& val = res.value();
      BOOST_CHECK_NE(val.lastMeasurementIndex, SIZE_MAX);
      BOOST_REQUIRE(val.fittedParameters);
      parameters = val.fittedParameters->parameters();
      BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
      // check the output status flags
      BOOST_CHECK(val.smoothed == expected_smoothed);
      BOOST_CHECK(val.reversed == expected_reversed);
      BOOST_CHECK(val.finished);
      BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
    }
    // fit w/ all hits in random order
    {
      auto shuffledSourceLinks = sourceLinks;
      std::shuffle(shuffledSourceLinks.begin(), shuffledSourceLinks.end(), rng);
      auto res = fitter.fit(shuffledSourceLinks.begin(),
                            shuffledSourceLinks.end(), start, options);
      BOOST_REQUIRE(res.ok());

      const auto& val = res.value();
      BOOST_CHECK_NE(val.lastMeasurementIndex, SIZE_MAX);
      BOOST_REQUIRE(val.fittedParameters);
      // check consistency w/ un-shuffled measurements
      CHECK_CLOSE_ABS(val.fittedParameters->parameters(), parameters, 1e-5);
      BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
      // check the output status flags
      BOOST_CHECK(val.smoothed == expected_smoothed);
      BOOST_CHECK(val.reversed == expected_reversed);
      BOOST_CHECK(val.finished);
      BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
    }
  }

  template <typename fitter_t, typename fitter_options_t, typename parameters_t>
  void test_ZeroFieldWithHole(const fitter_t& fitter, fitter_options_t options,
                              const parameters_t& start, Rng& rng,
                              const bool expected_reversed,
                              const bool expected_smoothed) const {
    auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                           resolutions, rng);
    const auto& sourceLinks = measurements.sourceLinks;
    BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

    // always keep the first and last measurement. leaving those in seems to not
    // count the respective surfaces as holes.
    for (size_t i = 1u; (i + 1u) < sourceLinks.size(); ++i) {
      // remove the i-th measurement
      auto withHole = sourceLinks;
      withHole.erase(std::next(withHole.begin(), i));
      BOOST_REQUIRE_EQUAL(withHole.size() + 1u, sourceLinks.size());
      BOOST_TEST_INFO("Removed measurement " << i);

      auto res = fitter.fit(withHole.begin(), withHole.end(), start, options);
      BOOST_REQUIRE(res.ok());

      const auto& val = res.value();
      BOOST_CHECK_NE(val.lastMeasurementIndex, SIZE_MAX);
      BOOST_CHECK(not val.fittedParameters);
      BOOST_CHECK_EQUAL(val.measurementStates, withHole.size());
      // check the output status flags
      BOOST_CHECK(val.smoothed == expected_smoothed);
      BOOST_CHECK(val.reversed == expected_reversed);
      BOOST_CHECK(val.finished);
      BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 1u);
    }
  }

  template <typename fitter_t, typename fitter_options_t, typename parameters_t>
  void test_ZeroFieldWithOutliers(const fitter_t& fitter,
                                  fitter_options_t options,
                                  const parameters_t& start, Rng& rng,
                                  const bool expected_reversed,
                                  const bool expected_smoothed) const {
    auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                           resolutions, rng);
    const auto& sourceLinks = measurements.sourceLinks;
    const auto& outlierSourceLinks = measurements.outlierSourceLinks;
    BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);
    BOOST_REQUIRE_EQUAL(outlierSourceLinks.size(), nMeasurements);

    for (size_t i = 0; i < sourceLinks.size(); ++i) {
      // replace the i-th measurement with an outlier
      auto withOutlier = sourceLinks;
      withOutlier[i] = outlierSourceLinks[i];
      BOOST_REQUIRE_EQUAL(withOutlier.size(), sourceLinks.size());
      BOOST_TEST_INFO("Replaced measurement " << i << " with outlier");

      auto res =
          fitter.fit(withOutlier.begin(), withOutlier.end(), start, options);
      BOOST_REQUIRE(res.ok());

      const auto& val = res.value();
      BOOST_CHECK_NE(val.lastMeasurementIndex, SIZE_MAX);
      // count the number of outliers
      size_t nOutliers = 0;
      val.fittedStates.visitBackwards(
          val.lastMeasurementIndex, [&nOutliers](const auto& state) {
            nOutliers +=
                state.typeFlags().test(Acts::TrackStateFlag::OutlierFlag);
          });
      BOOST_CHECK_EQUAL(nOutliers, 1u);
      BOOST_CHECK(not val.fittedParameters);
      BOOST_CHECK_EQUAL(val.measurementStates, withOutlier.size() - 1u);
      // check the output status flags
      BOOST_CHECK(val.smoothed == expected_smoothed);
      BOOST_CHECK(val.reversed == expected_reversed);
      BOOST_CHECK(val.finished);
      BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
    }
  }

  template <typename fitter_t, typename fitter_options_t, typename parameters_t>
  void test_ZeroFieldWithReverseFiltering(const fitter_t& fitter,
                                          fitter_options_t options,
                                          const parameters_t& start, Rng& rng,
                                          const bool expected_reversed,
                                          const bool expected_smoothed) const {
    auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                           resolutions, rng);
    const auto& sourceLinks = measurements.sourceLinks;
    const auto& outlierSourceLinks = measurements.outlierSourceLinks;
    BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);
    BOOST_REQUIRE_EQUAL(outlierSourceLinks.size(), nMeasurements);
    // create a boundless target surface near the tracker exit
    Acts::Vector3 center(3._m, 0., 0.);
    Acts::Vector3 normal(1., 0., 0.);
    auto targetSurface =
        Acts::Surface::makeShared<Acts::PlaneSurface>(center, normal);

    options.referenceSurface = targetSurface.get();

    auto res =
        fitter.fit(sourceLinks.begin(), sourceLinks.end(), start, options);
    BOOST_REQUIRE(res.ok());
    const auto& val = res.value();

    // Track of 1 GeV with a threshold set at 0.1 GeV, reversed filtering should
    // not be used
    BOOST_CHECK(val.smoothed == expected_smoothed);
    BOOST_CHECK(val.reversed == expected_reversed);
    BOOST_CHECK(val.finished);
  }

  // TODO this is not really Kalman fitter specific. is probably better tested
  // with a synthetic trajectory.
  template <typename fitter_t, typename fitter_options_t, typename parameters_t>
  void test_GlobalCovariance(const fitter_t& fitter, fitter_options_t options,
                             const parameters_t& start, Rng& rng) const {
    auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                           resolutions, rng);
    const auto& sourceLinks = measurements.sourceLinks;
    BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

    auto res =
        fitter.fit(sourceLinks.begin(), sourceLinks.end(), start, options);
    BOOST_REQUIRE(res.ok());

    // Calculate global track parameters covariance matrix
    const auto& val = res.value();
    auto [trackParamsCov, stateRowIndices] =
        Acts::detail::globalTrackParametersCovariance(val.fittedStates,
                                                      val.lastMeasurementIndex);
    BOOST_CHECK_EQUAL(trackParamsCov.rows(),
                      sourceLinks.size() * Acts::eBoundSize);
    BOOST_CHECK_EQUAL(stateRowIndices.size(), sourceLinks.size());
    // Each smoothed track state will have eBoundSize rows/cols in the global
    // covariance. stateRowIndices is a map of the starting row/index with the
    // state tip as the key. Thus, the last track state (i.e. the state
    // corresponding val.lastMeasurementIndex) has a starting row/index =
    // eBoundSize * (nMeasurements - 1), i.e. 6*(6-1) = 30.
    BOOST_CHECK_EQUAL(stateRowIndices.at(val.lastMeasurementIndex),
                      Acts::eBoundSize * (nMeasurements - 1));
  }
};
