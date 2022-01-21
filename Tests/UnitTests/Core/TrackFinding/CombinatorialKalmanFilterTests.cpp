// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <random>
#include <vector>

namespace {

using namespace Acts::Test;
using namespace Acts::UnitLiterals;

struct Detector {
  // expected number of measurements for the given detector
  size_t numMeasurements = 6u;

  // geometry
  CubicTrackingGeometry store;
  std::shared_ptr<const Acts::TrackingGeometry> geometry;

  // resolutions
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

  Detector(const Acts::GeometryContext geoCtx)
      : store(geoCtx), geometry(store()) {}
};

/// The map(-like) container accessor
template <typename container_t>
struct TestContainerAccessor {
  using Container = container_t;
  using Key = typename container_t::key_type;
  using Value = typename container_t::mapped_type;
  using Iterator = typename container_t::const_iterator;

  // pointer to the container
  const Container* container = nullptr;

  // count the number of elements with requested key
  size_t count(const Key& key) const {
    assert(container != nullptr);
    return container->count(key);
  }

  // get the range of elements with requested key
  std::pair<Iterator, Iterator> range(const Key& key) const {
    assert(container != nullptr);
    return container->equal_range(key);
  }

  // get the element using the iterator
  const Value& at(const Iterator& it) const { return (*it).second; }
};

struct Fixture {
  using StraightPropagator =
      Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
  using ConstantFieldStepper = Acts::EigenStepper<>;
  using ConstantFieldPropagator =
      Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

  using KalmanUpdater = Acts::GainMatrixUpdater;
  using KalmanSmoother = Acts::GainMatrixSmoother;
  using CombinatorialKalmanFilter =
      Acts::CombinatorialKalmanFilter<ConstantFieldPropagator>;
  using TestSourceLinkContainer =
      std::unordered_multimap<Acts::GeometryIdentifier, TestSourceLink>;
  using TestSourceLinkAccessor = TestContainerAccessor<TestSourceLinkContainer>;
  using CombinatorialKalmanFilterOptions =
      Acts::CombinatorialKalmanFilterOptions<TestSourceLinkAccessor>;

  KalmanUpdater kfUpdater;
  KalmanSmoother kfSmoother;

  Acts::GeometryContext geoCtx;
  Acts::MagneticFieldContext magCtx;
  Acts::CalibrationContext calCtx;

  Detector detector;

  // track parameters before and after the detector
  std::vector<Acts::CurvilinearTrackParameters> startParameters;
  std::vector<Acts::CurvilinearTrackParameters> endParameters;

  // generated measurements
  TestSourceLinkContainer sourceLinks;

  // CKF implementation to be tested
  CombinatorialKalmanFilter ckf;
  // configuration for the measurement selector
  Acts::MeasurementSelector::Config measurementSelectorCfg = {
      // global default: no chi2 cut, only one measurement per surface
      {Acts::GeometryIdentifier(),
       {{}, {std::numeric_limits<double>::max()}, {1u}}},
  };

  Acts::MeasurementSelector measSel{measurementSelectorCfg};

  Acts::CombinatorialKalmanFilterExtensions getExtensions() const {
    Acts::CombinatorialKalmanFilterExtensions extensions;
    extensions.calibrator.connect<&testSourceLinkCalibrator>();
    extensions.updater.connect<&KalmanUpdater::operator()>(&kfUpdater);
    extensions.smoother.connect<&KalmanSmoother::operator()>(&kfSmoother);
    extensions.measurementSelector.connect<&Acts::MeasurementSelector::select>(
        &measSel);
    return extensions;
  }

  std::unique_ptr<const Acts::Logger> logger;

  Fixture(double bz)
      : detector(geoCtx),
        ckf(makeConstantFieldPropagator(detector.geometry, bz)),
        logger(Acts::getDefaultLogger("CkfTest", Acts::Logging::INFO)) {
    // construct initial parameters
    // create common covariance matrix from reasonable standard deviations
    Acts::BoundVector stddev;
    stddev[Acts::eBoundLoc0] = 100_um;
    stddev[Acts::eBoundLoc1] = 100_um;
    stddev[Acts::eBoundTime] = 25_ns;
    stddev[Acts::eBoundPhi] = 2_degree;
    stddev[Acts::eBoundTheta] = 2_degree;
    stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
    Acts::BoundSymMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
    // all tracks close to the transverse plane along the x axis w/ small
    // variations in position, direction.
    Acts::Vector4 mStartPos0(-3_m, 0.0, 0.0, 1_ns);
    Acts::Vector4 mStartPos1(-3_m, -15_mm, -15_mm, 2_ns);
    Acts::Vector4 mStartPos2(-3_m, 15_mm, 15_mm, -1_ns);
    startParameters = {
        {mStartPos0, 0_degree, 90_degree, 1_GeV, 1_e, cov},
        {mStartPos1, -1_degree, 91_degree, 1_GeV, 1_e, cov},
        {mStartPos2, 1_degree, 89_degree, 1_GeV, -1_e, cov},
    };
    Acts::Vector4 mEndPos0(3_m, 0.0, 0.0, 1_ns);
    Acts::Vector4 mEndPos1(3_m, -100_mm, -100_mm, 2_ns);
    Acts::Vector4 mEndPos2(3_m, 100_mm, 100_mm, -1_ns);
    endParameters = {
        {mEndPos0, 0_degree, 90_degree, 1_GeV, 1_e, cov * 100},
        {mEndPos1, -1_degree, 91_degree, 1_GeV, 1_e, cov * 100},
        {mEndPos2, 1_degree, 89_degree, 1_GeV, -1_e, cov * 100},
    };

    // create some measurements
    auto measPropagator = makeStraightPropagator(detector.geometry);
    std::default_random_engine rng(421235);
    for (size_t trackId = 0u; trackId < startParameters.size(); ++trackId) {
      auto measurements = createMeasurements(
          measPropagator, geoCtx, magCtx, startParameters[trackId],
          detector.resolutions, rng, trackId);
      for (auto& sl : measurements.sourceLinks) {
        sourceLinks.emplace(sl.geometryId(), std::move(sl));
      }
    }
  }

  // Construct a straight-line propagator.
  static StraightPropagator makeStraightPropagator(
      std::shared_ptr<const Acts::TrackingGeometry> geo) {
    Acts::Navigator::Config cfg{geo};
    cfg.resolvePassive = false;
    cfg.resolveMaterial = true;
    cfg.resolveSensitive = true;
    Acts::Navigator navigator{cfg};
    Acts::StraightLineStepper stepper;
    return StraightPropagator(std::move(stepper), std::move(navigator));
  }

  // Construct a propagator using a constant magnetic field along z.
  static ConstantFieldPropagator makeConstantFieldPropagator(
      std::shared_ptr<const Acts::TrackingGeometry> geo, double bz) {
    Acts::Navigator::Config cfg{geo};
    cfg.resolvePassive = false;
    cfg.resolveMaterial = true;
    cfg.resolveSensitive = true;
    Acts::Navigator navigator{cfg};
    auto field =
        std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, bz));
    ConstantFieldStepper stepper(std::move(field));
    return ConstantFieldPropagator(std::move(stepper), std::move(navigator));
  }

  CombinatorialKalmanFilterOptions makeCkfOptions() const {
    return CombinatorialKalmanFilterOptions(
        geoCtx, magCtx, calCtx, TestSourceLinkAccessor(), getExtensions(),
        Acts::LoggerWrapper{*logger}, Acts::PropagatorPlainOptions());
  }
};

}  // namespace

BOOST_AUTO_TEST_SUITE(TrackFindingCombinatorialKalmanFilter)

BOOST_AUTO_TEST_CASE(ZeroFieldForward) {
  Fixture f(0_T);

  auto options = f.makeCkfOptions();
  // this is the default option. set anyways for consistency
  options.propagatorPlainOptions.direction = Acts::forward;
  // Construct a plane surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Vector3{-3_m, 0., 0.}, Acts::Vector3{1., 0., 0});
  // Set the target surface
  options.referenceSurface = &(*pSurface);

  // run the CKF for all initial track states
  auto results = f.ckf.findTracks(f.sourceLinks, f.startParameters, options);
  // There should be three track finding results with three initial track states
  BOOST_CHECK_EQUAL(results.size(), 3u);

  // check the found tracks
  for (size_t trackId = 0u; trackId < f.startParameters.size(); ++trackId) {
    const auto& params = f.startParameters[trackId];
    BOOST_TEST_INFO("initial parameters before detector:\n" << params);

    auto& res = results[trackId];
    if (not res.ok()) {
      BOOST_TEST_INFO(res.error() << " " << res.error().message());
    }
    BOOST_REQUIRE(res.ok());

    auto val = *res;
    BOOST_CHECK_EQUAL(val.fittedStates.size(), f.detector.numMeasurements);

    // with the given measurement selection cuts, only one trajectory for the
    // given input parameters should be found.
    BOOST_CHECK_EQUAL(val.lastMeasurementIndices.size(), 1u);
    // check purity of first found track
    // find the number of hits not originating from the right track
    size_t numHits = 0u;
    size_t nummismatchedHits = 0u;
    val.fittedStates.visitBackwards(
        val.lastMeasurementIndices.front(), [&](const auto& trackState) {
          numHits += 1u;
          const auto& sl =
              static_cast<const TestSourceLink&>(trackState.uncalibrated());
          nummismatchedHits += (trackId != sl.sourceId);
        });

    BOOST_CHECK_EQUAL(numHits, f.detector.numMeasurements);
    BOOST_CHECK_EQUAL(nummismatchedHits, 0u);
  }
}

BOOST_AUTO_TEST_CASE(ZeroFieldBackward) {
  Fixture f(0_T);

  auto options = f.makeCkfOptions();
  options.propagatorPlainOptions.direction = Acts::backward;
  // Construct a plane surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Vector3{3_m, 0., 0.}, Acts::Vector3{1., 0., 0});
  // Set the target surface
  options.referenceSurface = &(*pSurface);

  // run the CKF for all initial track states
  auto results = f.ckf.findTracks(f.sourceLinks, f.endParameters, options);
  // There should be three found tracks with three initial track states
  BOOST_CHECK_EQUAL(results.size(), 3u);

  // check the found tracks
  for (size_t trackId = 0u; trackId < f.endParameters.size(); ++trackId) {
    const auto& params = f.endParameters[trackId];
    BOOST_TEST_INFO("initial parameters after detector:\n" << params);

    auto& res = results[trackId];
    if (not res.ok()) {
      BOOST_TEST_INFO(res.error() << " " << res.error().message());
    }
    BOOST_REQUIRE(res.ok());

    auto val = *res;
    BOOST_CHECK_EQUAL(val.fittedStates.size(), f.detector.numMeasurements);
    // with the given measurement selection cuts, only one trajectory for the
    // given input parameters should be found.
    BOOST_CHECK_EQUAL(val.lastMeasurementIndices.size(), 1u);
    // check purity of first found track
    // find the number of hits not originating from the right track
    size_t numHits = 0u;
    size_t nummismatchedHits = 0u;
    val.fittedStates.visitBackwards(
        val.lastMeasurementIndices.front(), [&](const auto& trackState) {
          numHits += 1u;
          nummismatchedHits += (trackId != static_cast<const TestSourceLink&>(
                                               trackState.uncalibrated())
                                               .sourceId);
        });
    BOOST_CHECK_EQUAL(numHits, f.detector.numMeasurements);
    BOOST_CHECK_EQUAL(nummismatchedHits, 0u);
  }
}

BOOST_AUTO_TEST_SUITE_END()
