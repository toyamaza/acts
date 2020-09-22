// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <iostream>
#include <numeric>
#include <random>

using std::cout;
using std::endl;

namespace Acts {
namespace Test {

GeometryContext gctx;

using SourceLink = MinimalSourceLink;
using Parameters = BoundVector;
using Covariance = BoundSymMatrix;

CurvilinearTrackParameters make_params() {
  // generate arbitrary positive, definite matrix
  Covariance rnd = Covariance::Random();
  Covariance cov = rnd.transpose() * rnd;
  return {Vector4D(0, 0, 1, 0), Vector3D(1, 10, 40), 1000, -1, cov};
}

using ParVec_t = BoundTrackParameters::ParametersVector;
using CovMat_t = BoundTrackParameters::CovarianceMatrix;

struct TestTrackState {
  SourceLink sourceLink;
  std::optional<Measurement<SourceLink, BoundIndices, eBoundLoc0, eBoundLoc1,
                            eBoundQOverP>>
      meas3d;
  std::optional<Measurement<SourceLink, BoundIndices, eBoundLoc0, eBoundLoc1>>
      meas2d;
  std::optional<BoundTrackParameters> predicted;
  std::optional<BoundTrackParameters> filtered;
  std::optional<BoundTrackParameters> smoothed;
  CovMat_t jacobian;
  double chi2;
  double pathLength;
};

/// @brief Fills a @c TrackStateProxy object
///
/// @tparam track_state_t Type of the TrackStateProxy
///
/// @param [in, out] ts TrackStateProxy which is filled
/// @param [in] mask Specifies which components are filled
/// @param [in] dim Dimension of the measurement
///
/// @return Tuple containing a @c TestTrackState and the @c FittableMeasurement
/// that were generated in this function
template <typename track_state_t>
auto fillTrackState(track_state_t& ts, TrackStatePropMask mask,
                    size_t dim = 3) {
  auto plane = Surface::makeShared<PlaneSurface>(Vector3D{0., 0., 0.},
                                                 Vector3D{0., 0., 1.});

  std::unique_ptr<FittableMeasurement<SourceLink>> fm;
  TestTrackState pc;

  if (dim == 3) {
    ActsMatrixD<3, 3> mCov;
    mCov.setRandom();

    Vector3D mPar;
    mPar.setRandom();
    Measurement<SourceLink, BoundIndices, eBoundLoc0, eBoundLoc1, eBoundQOverP>
        meas{plane, {}, mCov, mPar[0], mPar[1], mPar[2]};

    fm = std::make_unique<FittableMeasurement<SourceLink>>(meas);

    SourceLink sourceLink{fm.get()};
    pc.sourceLink = sourceLink;
    if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Uncalibrated)) {
      ts.uncalibrated() = sourceLink;
    }

    // "calibrate", keep original source link (stack address)
    pc.meas3d = {meas.referenceObject().getSharedPtr(),
                 sourceLink,
                 meas.covariance(),
                 meas.parameters()[0],
                 meas.parameters()[1],
                 meas.parameters()[2]};
    if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Calibrated)) {
      ts.setCalibrated(*pc.meas3d);
    }
  } else if (dim == 2) {
    ActsMatrixD<2, 2> mCov;
    mCov.setRandom();

    Vector2D mPar;
    mPar.setRandom();
    Measurement<SourceLink, BoundIndices, eBoundLoc0, eBoundLoc1> meas{
        plane, {}, mCov, mPar[0], mPar[1]};

    fm = std::make_unique<FittableMeasurement<SourceLink>>(meas);

    SourceLink sourceLink{fm.get()};
    pc.sourceLink = sourceLink;
    if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Uncalibrated)) {
      ts.uncalibrated() = sourceLink;
    }

    // "calibrate", keep original source link (stack address)
    pc.meas2d = {meas.referenceObject().getSharedPtr(), sourceLink,
                 meas.covariance(), meas.parameters()[0], meas.parameters()[1]};
    if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Calibrated)) {
      ts.setCalibrated(*pc.meas2d);
    }
  } else {
    throw std::runtime_error("wrong dim");
  }

  // add parameters

  // predicted
  ParVec_t predPar;
  predPar << 1, 2, M_PI / 4., M_PI / 2., 5, 0.;
  predPar.template head<2>().setRandom();

  CovMat_t predCov;
  predCov.setRandom();

  BoundTrackParameters pred(plane, predPar, predCov);
  pc.predicted = pred;
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Predicted)) {
    ts.predicted() = pred.parameters();
    ts.predictedCovariance() = *pred.covariance();
  }

  // filtered
  ParVec_t filtPar;
  filtPar << 6, 7, M_PI / 4., M_PI / 2., 10, 0.;
  filtPar.template head<2>().setRandom();

  CovMat_t filtCov;
  filtCov.setRandom();

  BoundTrackParameters filt(plane, filtPar, filtCov);
  pc.filtered = filt;
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Filtered)) {
    ts.filtered() = filt.parameters();
    ts.filteredCovariance() = *filt.covariance();
  }

  // smoothed
  ParVec_t smotPar;
  smotPar << 11, 12, M_PI / 4., M_PI / 2., 15, 0.;
  smotPar.template head<2>().setRandom();

  CovMat_t smotCov;
  smotCov.setRandom();

  BoundTrackParameters smot(plane, smotPar, smotCov);
  pc.smoothed = smot;
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Smoothed)) {
    ts.smoothed() = smot.parameters();
    ts.smoothedCovariance() = *smot.covariance();
  }

  // make jacobian
  CovMat_t jac;
  jac.setRandom();
  pc.jacobian = jac;
  if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Jacobian)) {
    ts.jacobian() = jac;
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(1.0, 100.0);
  pc.chi2 = dis(gen);
  pc.pathLength = dis(gen);
  ts.chi2() = pc.chi2;
  ts.pathLength() = pc.pathLength;

  return std::make_tuple(pc, std::move(fm));
}

BOOST_AUTO_TEST_CASE(multitrajectory_build) {
  MultiTrajectory<SourceLink> t;
  TrackStatePropMask mask = TrackStatePropMask::Predicted;

  // construct trajectory w/ multiple components
  auto i0 = t.addTrackState(mask);
  // trajectory bifurcates here into multiple hypotheses
  auto i1a = t.addTrackState(mask, i0);
  auto i1b = t.addTrackState(mask, i0);
  auto i2a = t.addTrackState(mask, i1a);
  auto i2b = t.addTrackState(mask, i1b);

  // print each trajectory component
  std::vector<size_t> act;
  auto collect = [&](auto p) {
    act.push_back(p.index());
    BOOST_CHECK(!p.hasUncalibrated());
    BOOST_CHECK(!p.hasCalibrated());
    BOOST_CHECK(!p.hasFiltered());
    BOOST_CHECK(!p.hasSmoothed());
    BOOST_CHECK(!p.hasJacobian());
    BOOST_CHECK(!p.hasProjector());
  };

  std::vector<size_t> exp = {i2a, i1a, i0};
  t.visitBackwards(i2a, collect);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());

  act.clear();
  exp = {i2b, i1b, i0};
  t.visitBackwards(i2b, collect);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());

  act.clear();
  t.applyBackwards(i2b, collect);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());
}

BOOST_AUTO_TEST_CASE(visit_apply_abort) {
  MultiTrajectory<SourceLink> t;
  TrackStatePropMask mask = TrackStatePropMask::Predicted;

  // construct trajectory with three components
  auto i0 = t.addTrackState(mask);
  auto i1 = t.addTrackState(mask, i0);
  auto i2 = t.addTrackState(mask, i1);

  size_t n = 0;
  t.applyBackwards(i2, [&](const auto&) {
    n++;
    return false;
  });
  BOOST_CHECK_EQUAL(n, 1u);

  n = 0;
  t.applyBackwards(i2, [&](const auto& ts) {
    n++;
    if (ts.index() == i1) {
      return false;
    }
    return true;
  });
  BOOST_CHECK_EQUAL(n, 2u);

  n = 0;
  t.applyBackwards(i2, [&](const auto&) {
    n++;
    return true;
  });
  BOOST_CHECK_EQUAL(n, 3u);
}

BOOST_AUTO_TEST_CASE(trackstate_add_bitmask_operators) {
  using PM = TrackStatePropMask;
  auto bs1 = PM::Uncalibrated;

  BOOST_CHECK(ACTS_CHECK_BIT(bs1, PM::Uncalibrated));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs1, PM::Calibrated));

  auto bs2 = PM::Calibrated;

  BOOST_CHECK(!ACTS_CHECK_BIT(bs2, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(bs2, PM::Calibrated));

  auto bs3 = PM::Calibrated | PM::Uncalibrated;

  BOOST_CHECK(ACTS_CHECK_BIT(bs3, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(bs3, PM::Calibrated));

  BOOST_CHECK(ACTS_CHECK_BIT(PM::All, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(PM::All, PM::Calibrated));

  auto bs4 = PM::Predicted | PM::Jacobian | PM::Uncalibrated;
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Predicted));
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Uncalibrated));
  BOOST_CHECK(ACTS_CHECK_BIT(bs4, PM::Jacobian));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs4, PM::Calibrated));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs4, PM::Filtered));
  BOOST_CHECK(!ACTS_CHECK_BIT(bs4, PM::Smoothed));

  auto cnv = [](auto a) -> std::bitset<8> {
    return static_cast<std::underlying_type<PM>::type>(a);
  };

  BOOST_CHECK(cnv(PM::All).all());    // all ones
  BOOST_CHECK(cnv(PM::None).none());  // all zeros

  // test orthogonality
  std::array<PM, 6> values{PM::Predicted, PM::Filtered,     PM::Smoothed,
                           PM::Jacobian,  PM::Uncalibrated, PM::Calibrated};
  for (size_t i = 0; i < values.size(); i++) {
    for (size_t j = 0; j < values.size(); j++) {
      PM a = values[i];
      PM b = values[j];

      if (i == j) {
        BOOST_CHECK(cnv(a & b).count() == 1);
      } else {
        BOOST_CHECK(cnv(a & b).none());
      }
    }
  }

  BOOST_CHECK(cnv(PM::Predicted ^ PM::Filtered).count() == 2);
  BOOST_CHECK(cnv(PM::Predicted ^ PM::Predicted).none());
  BOOST_CHECK(~(PM::Predicted | PM::Calibrated) ==
              (PM::All ^ PM::Predicted ^ PM::Calibrated));

  PM base = PM::None;
  BOOST_CHECK(cnv(base) == 0);

  base &= PM::Filtered;
  BOOST_CHECK(cnv(base) == 0);

  base |= PM::Filtered;
  BOOST_CHECK(base == PM::Filtered);

  base |= PM::Calibrated;
  BOOST_CHECK(base == (PM::Filtered | PM::Calibrated));

  base ^= PM::All;
  BOOST_CHECK(base == ~(PM::Filtered | PM::Calibrated));
}

BOOST_AUTO_TEST_CASE(trackstate_add_bitmask_method) {
  using PM = TrackStatePropMask;
  MultiTrajectory<SourceLink> t;

  auto ts = t.getTrackState(t.addTrackState(PM::All));
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK(ts.hasSmoothed());
  BOOST_CHECK(ts.hasUncalibrated());
  BOOST_CHECK(ts.hasCalibrated());
  BOOST_CHECK(ts.hasProjector());
  BOOST_CHECK(ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::None));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Predicted));
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Filtered));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Smoothed));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Uncalibrated));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Calibrated));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(ts.hasCalibrated());
  BOOST_CHECK(ts.hasProjector());
  BOOST_CHECK(!ts.hasJacobian());

  ts = t.getTrackState(t.addTrackState(PM::Jacobian));
  BOOST_CHECK(!ts.hasPredicted());
  BOOST_CHECK(!ts.hasFiltered());
  BOOST_CHECK(!ts.hasSmoothed());
  BOOST_CHECK(!ts.hasUncalibrated());
  BOOST_CHECK(!ts.hasCalibrated());
  BOOST_CHECK(!ts.hasProjector());
  BOOST_CHECK(ts.hasJacobian());
}

BOOST_AUTO_TEST_CASE(trackstate_proxy_cross_talk) {
  // assert expected "cross-talk" between trackstate proxies

  MultiTrajectory<SourceLink> t;
  size_t index = t.addTrackState();
  auto tso = t.getTrackState(index);
  auto [pc, fm] = fillTrackState(tso, TrackStatePropMask::All);

  const auto& ct = t;
  auto cts = ct.getTrackState(0);
  auto ts = t.getTrackState(0);

  // assert expected value of chi2 and path length
  BOOST_CHECK_EQUAL(cts.chi2(), pc.chi2);
  BOOST_CHECK_EQUAL(ts.chi2(), pc.chi2);
  BOOST_CHECK_EQUAL(cts.pathLength(), pc.pathLength);
  BOOST_CHECK_EQUAL(ts.pathLength(), pc.pathLength);

  ParVec_t v;
  CovMat_t cov;

  v.setRandom();
  ts.predicted() = v;
  BOOST_CHECK_EQUAL(cts.predicted(), v);
  cov.setRandom();
  ts.predictedCovariance() = cov;
  BOOST_CHECK_EQUAL(cts.predictedCovariance(), cov);

  v.setRandom();
  ts.filtered() = v;
  BOOST_CHECK_EQUAL(cts.filtered(), v);
  cov.setRandom();
  ts.filteredCovariance() = cov;
  BOOST_CHECK_EQUAL(cts.filteredCovariance(), cov);

  v.setRandom();
  ts.smoothed() = v;
  BOOST_CHECK_EQUAL(cts.smoothed(), v);
  cov.setRandom();
  ts.smoothedCovariance() = cov;
  BOOST_CHECK_EQUAL(cts.smoothedCovariance(), cov);

  // make copy of fm
  auto fm2 = std::make_unique<FittableMeasurement<SourceLink>>(*fm);
  SourceLink sourceLink2{fm2.get()};
  ts.uncalibrated() = sourceLink2;
  BOOST_CHECK_EQUAL(cts.uncalibrated(), sourceLink2);
  BOOST_CHECK_NE(cts.uncalibrated(), SourceLink{fm.get()});

  CovMat_t newMeasCov;
  newMeasCov.setRandom();
  ts.calibratedCovariance() = newMeasCov;
  BOOST_CHECK_EQUAL(cts.calibratedCovariance(), newMeasCov);

  ParVec_t newMeasPar;
  newMeasPar.setRandom();
  ts.calibrated() = newMeasPar;
  BOOST_CHECK_EQUAL(cts.calibrated(), newMeasPar);

  size_t measdim = ts.effectiveCalibrated().rows();

  ActsMatrixXd eff{measdim, measdim};
  eff.setRandom();
  ts.effectiveCalibratedCovariance() = eff;
  BOOST_CHECK_EQUAL(cts.effectiveCalibratedCovariance(), eff);
  newMeasCov.topLeftCorner(eff.rows(), eff.rows()) = eff;
  BOOST_CHECK_EQUAL(cts.calibratedCovariance(), newMeasCov);

  CovMat_t jac;
  jac.setRandom();
  ts.jacobian() = jac;
  BOOST_CHECK_EQUAL(cts.jacobian(), jac);

  ts.chi2() = 98;
  BOOST_CHECK_EQUAL(cts.chi2(), 98u);

  ts.pathLength() = 66;
  BOOST_CHECK_EQUAL(cts.pathLength(), 66u);
}

BOOST_AUTO_TEST_CASE(trackstate_reassignment) {
  constexpr size_t maxmeasdim = MultiTrajectory<SourceLink>::MeasurementSizeMax;

  MultiTrajectory<SourceLink> t;
  size_t index = t.addTrackState();
  auto tso = t.getTrackState(index);
  auto [pc, fm] = fillTrackState(tso, TrackStatePropMask::All);

  auto ts = t.getTrackState(0);

  // assert measdim and contents of original measurement (just to be safe)
  BOOST_CHECK_EQUAL(ts.calibratedSize(), pc.meas3d->size());
  BOOST_CHECK_EQUAL(ts.effectiveCalibrated(), pc.meas3d->parameters());
  BOOST_CHECK_EQUAL(ts.effectiveCalibratedCovariance(),
                    pc.meas3d->covariance());
  BOOST_CHECK_EQUAL(ts.effectiveProjector(), pc.meas3d->projector());

  // create new measurement
  SymMatrix2D mCov;
  mCov.setRandom();
  Vector2D mPar;
  mPar.setRandom();
  Measurement<SourceLink, BoundIndices, eBoundLoc0, eBoundLoc1> m2{
      pc.meas3d->referenceObject().getSharedPtr(), {}, mCov, mPar[0], mPar[1]};

  ts.setCalibrated(m2);

  BOOST_CHECK_EQUAL(ts.calibratedSize(), 2u);
  BOOST_CHECK_EQUAL(ts.effectiveCalibrated(), mPar);
  BOOST_CHECK_EQUAL(ts.effectiveCalibratedCovariance(), mCov);
  BOOST_CHECK_EQUAL(ts.effectiveProjector(), m2.projector());

  // check if overallocated part is zeroed correctly
  ActsVectorD<maxmeasdim> mParFull;
  mParFull.setZero();
  mParFull.head(2) = mPar;
  BOOST_CHECK_EQUAL(ts.calibrated(), mParFull);

  BoundSymMatrix mCovFull;
  mCovFull.setZero();
  mCovFull.topLeftCorner(2, 2) = mCov;
  BOOST_CHECK_EQUAL(ts.calibratedCovariance(), mCovFull);

  ActsMatrixD<maxmeasdim, eBoundSize> projFull;
  projFull.setZero();
  projFull.topLeftCorner(m2.size(), eBoundSize) = m2.projector();
  BOOST_CHECK_EQUAL(ts.projector(), projFull);
}

BOOST_AUTO_TEST_CASE(storage_consistency) {
  MultiTrajectory<SourceLink> t;
  size_t index = t.addTrackState();
  auto ts = t.getTrackState(index);
  auto [pc, fm] = fillTrackState(ts, TrackStatePropMask::All);

  // now investigate the proxy
  // parameters
  BOOST_CHECK(ts.hasPredicted());
  BOOST_CHECK_EQUAL(pc.predicted->parameters(), ts.predicted());
  BOOST_CHECK_EQUAL(*pc.predicted->covariance(), ts.predictedCovariance());

  BOOST_CHECK(ts.hasFiltered());
  BOOST_CHECK_EQUAL(pc.filtered->parameters(), ts.filtered());
  BOOST_CHECK_EQUAL(*pc.filtered->covariance(), ts.filteredCovariance());

  BOOST_CHECK(ts.hasSmoothed());
  BOOST_CHECK_EQUAL(pc.smoothed->parameters(), ts.smoothed());
  BOOST_CHECK_EQUAL(*pc.smoothed->covariance(), ts.smoothedCovariance());

  BOOST_CHECK_EQUAL(&ts.referenceSurface(), &pc.sourceLink.referenceSurface());

  BOOST_CHECK(ts.hasJacobian());
  BOOST_CHECK_EQUAL(ts.jacobian(), pc.jacobian);

  BOOST_CHECK(ts.hasProjector());
  BOOST_CHECK_EQUAL(ts.effectiveProjector(), pc.meas3d->projector());
  // measurement properties
  BOOST_CHECK(ts.hasCalibrated());
  BOOST_CHECK_EQUAL(pc.meas3d->parameters(), ts.effectiveCalibrated());
  ParVec_t mParFull;
  mParFull.setZero();
  mParFull.head(pc.meas3d->size()) = pc.meas3d->parameters();
  BOOST_CHECK_EQUAL(mParFull, ts.calibrated());

  BOOST_CHECK_EQUAL(pc.meas3d->covariance(),
                    ts.effectiveCalibratedCovariance());
  CovMat_t mCovFull;
  mCovFull.setZero();
  mCovFull.topLeftCorner(pc.meas3d->size(), pc.meas3d->size()) =
      pc.meas3d->covariance();
  BOOST_CHECK_EQUAL(mCovFull, ts.calibratedCovariance());

  // calibrated links to original measurement
  BOOST_CHECK_EQUAL(pc.meas3d->sourceLink(), ts.calibratedSourceLink());

  // uncalibrated **is** a SourceLink
  BOOST_CHECK(ts.hasUncalibrated());
  BOOST_CHECK_EQUAL(pc.meas3d->sourceLink(), ts.uncalibrated());

  // full projector, should be exactly equal
  ActsMatrixD<MultiTrajectory<SourceLink>::MeasurementSizeMax, eBoundSize>
      fullProj;
  fullProj.setZero();
  fullProj.topLeftCorner(pc.meas3d->size(), eBoundSize) =
      pc.meas3d->projector();
  BOOST_CHECK_EQUAL(ts.projector(), fullProj);

  // projector with dynamic rows
  // should be exactly equal
  BOOST_CHECK_EQUAL(ts.effectiveProjector(), pc.meas3d->projector());
}

BOOST_AUTO_TEST_CASE(add_trackstate_allocations) {
  MultiTrajectory<SourceLink> t;

  // this should allocate for all the components in the trackstate, plus
  // filtered
  size_t i = t.addTrackState(
      TrackStatePropMask::Predicted | TrackStatePropMask::Filtered |
      TrackStatePropMask::Uncalibrated | TrackStatePropMask::Jacobian);
  auto tso = t.getTrackState(i);
  fillTrackState(tso, TrackStatePropMask::Predicted);
  fillTrackState(tso, TrackStatePropMask::Filtered);
  fillTrackState(tso, TrackStatePropMask::Uncalibrated);
  fillTrackState(tso, TrackStatePropMask::Jacobian);

  BOOST_CHECK(tso.hasPredicted());
  BOOST_CHECK(tso.hasFiltered());
  BOOST_CHECK(!tso.hasSmoothed());
  BOOST_CHECK(tso.hasUncalibrated());
  BOOST_CHECK(!tso.hasCalibrated());
  BOOST_CHECK(tso.hasJacobian());

  // remove some parts
}

BOOST_AUTO_TEST_CASE(trackstateproxy_getmask) {
  using PM = TrackStatePropMask;
  MultiTrajectory<SourceLink> mj;

  std::array<PM, 6> values{PM::Predicted, PM::Filtered,     PM::Smoothed,
                           PM::Jacobian,  PM::Uncalibrated, PM::Calibrated};

  PM all = std::accumulate(values.begin(), values.end(), PM::None,
                           [](auto a, auto b) { return a | b; });

  auto ts = mj.getTrackState(mj.addTrackState(PM::All));
  BOOST_CHECK(ts.getMask() == all);

  ts = mj.getTrackState(mj.addTrackState(PM::Filtered | PM::Calibrated));
  BOOST_CHECK(ts.getMask() == (PM::Filtered | PM::Calibrated));

  ts = mj.getTrackState(
      mj.addTrackState(PM::Filtered | PM::Smoothed | PM::Predicted));
  BOOST_CHECK(ts.getMask() == (PM::Filtered | PM::Smoothed | PM::Predicted));

  for (PM mask : values) {
    ts = mj.getTrackState(mj.addTrackState(mask));
    BOOST_CHECK(ts.getMask() == mask);
  }
}

BOOST_AUTO_TEST_CASE(trackstateproxy_copy) {
  using PM = TrackStatePropMask;
  MultiTrajectory<SourceLink> mj;
  auto mkts = [&](PM mask) { return mj.getTrackState(mj.addTrackState(mask)); };

  std::array<PM, 6> values{PM::Predicted, PM::Filtered,     PM::Smoothed,
                           PM::Jacobian,  PM::Uncalibrated, PM::Calibrated};

  // orthogonal ones

  for (PM a : values) {
    for (PM b : values) {
      auto tsa = mkts(a);
      auto tsb = mkts(b);
      // doesn't work
      if (a != b) {
        BOOST_CHECK_THROW(tsa.copyFrom(tsb), std::runtime_error);
        BOOST_CHECK_THROW(tsb.copyFrom(tsa), std::runtime_error);
      } else {
        tsa.copyFrom(tsb);
        tsb.copyFrom(tsa);
      }
    }
  }

  auto ts1 = mkts(PM::Filtered | PM::Predicted);  // this has both
  ts1.filtered().setRandom();
  ts1.filteredCovariance().setRandom();
  ts1.predicted().setRandom();
  ts1.predictedCovariance().setRandom();

  // ((src XOR dst) & src) == 0
  auto ts2 = mkts(PM::Predicted);
  ts2.predicted().setRandom();
  ts2.predictedCovariance().setRandom();

  // they are different before
  BOOST_CHECK(ts1.predicted() != ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() != ts2.predictedCovariance());

  // ts1 -> ts2 fails
  BOOST_CHECK_THROW(ts2.copyFrom(ts1), std::runtime_error);
  BOOST_CHECK(ts1.predicted() != ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() != ts2.predictedCovariance());

  // ts2 -> ts1 is ok
  ts1.copyFrom(ts2);
  BOOST_CHECK(ts1.predicted() == ts2.predicted());
  BOOST_CHECK(ts1.predictedCovariance() == ts2.predictedCovariance());

  size_t i0 = mj.addTrackState();
  size_t i1 = mj.addTrackState();
  ts1 = mj.getTrackState(i0);
  ts2 = mj.getTrackState(i1);
  auto [rts1, fm1] = fillTrackState(ts1, TrackStatePropMask::All, 2);
  auto [rts2, fm2] = fillTrackState(ts2, TrackStatePropMask::All, 3);

  auto ots1 = mkts(PM::All);
  auto ots2 = mkts(PM::All);
  // make full copy for later. We prove full copy works right below
  ots1.copyFrom(ts1);
  ots2.copyFrom(ts2);

  BOOST_CHECK_NE(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_NE(ts1.predictedCovariance(), ts2.predictedCovariance());
  BOOST_CHECK_NE(ts1.filtered(), ts2.filtered());
  BOOST_CHECK_NE(ts1.filteredCovariance(), ts2.filteredCovariance());
  BOOST_CHECK_NE(ts1.smoothed(), ts2.smoothed());
  BOOST_CHECK_NE(ts1.smoothedCovariance(), ts2.smoothedCovariance());

  BOOST_CHECK_NE(ts1.uncalibrated(), ts2.uncalibrated());

  BOOST_CHECK_NE(ts1.calibratedSourceLink(), ts2.calibratedSourceLink());
  BOOST_CHECK_NE(ts1.calibrated(), ts2.calibrated());
  BOOST_CHECK_NE(ts1.calibratedCovariance(), ts2.calibratedCovariance());
  BOOST_CHECK_NE(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_NE(ts1.projector(), ts2.projector());

  BOOST_CHECK_NE(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_NE(ts1.chi2(), ts2.chi2());
  BOOST_CHECK_NE(ts1.pathLength(), ts2.pathLength());
  BOOST_CHECK_NE(&ts1.referenceSurface(), &ts2.referenceSurface());

  ts1.copyFrom(ts2);

  BOOST_CHECK_EQUAL(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_EQUAL(ts1.predictedCovariance(), ts2.predictedCovariance());
  BOOST_CHECK_EQUAL(ts1.filtered(), ts2.filtered());
  BOOST_CHECK_EQUAL(ts1.filteredCovariance(), ts2.filteredCovariance());
  BOOST_CHECK_EQUAL(ts1.smoothed(), ts2.smoothed());
  BOOST_CHECK_EQUAL(ts1.smoothedCovariance(), ts2.smoothedCovariance());

  BOOST_CHECK_EQUAL(ts1.uncalibrated(), ts2.uncalibrated());

  BOOST_CHECK_EQUAL(ts1.calibratedSourceLink(), ts2.calibratedSourceLink());
  BOOST_CHECK_EQUAL(ts1.calibrated(), ts2.calibrated());
  BOOST_CHECK_EQUAL(ts1.calibratedCovariance(), ts2.calibratedCovariance());
  BOOST_CHECK_EQUAL(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_EQUAL(ts1.projector(), ts2.projector());

  BOOST_CHECK_EQUAL(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_EQUAL(ts1.chi2(), ts2.chi2());
  BOOST_CHECK_EQUAL(ts1.pathLength(), ts2.pathLength());
  BOOST_CHECK_EQUAL(&ts1.referenceSurface(), &ts2.referenceSurface());

  // full copy proven to work. now let's do partial copy
  ts2 = mkts(PM::Predicted | PM::Jacobian | PM::Calibrated);
  ts2.copyFrom(ots2, PM::Predicted | PM::Jacobian |
                         PM::Calibrated);  // copy into empty ts, only copy some
  ts1.copyFrom(ots1);                      // reset to original
  // is different again
  BOOST_CHECK_NE(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_NE(ts1.predictedCovariance(), ts2.predictedCovariance());

  BOOST_CHECK_NE(ts1.calibratedSourceLink(), ts2.calibratedSourceLink());
  BOOST_CHECK_NE(ts1.calibrated(), ts2.calibrated());
  BOOST_CHECK_NE(ts1.calibratedCovariance(), ts2.calibratedCovariance());
  BOOST_CHECK_NE(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_NE(ts1.projector(), ts2.projector());

  BOOST_CHECK_NE(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_NE(ts1.chi2(), ts2.chi2());
  BOOST_CHECK_NE(ts1.pathLength(), ts2.pathLength());
  BOOST_CHECK_NE(&ts1.referenceSurface(), &ts2.referenceSurface());

  ts1.copyFrom(ts2);

  // some components are same now
  BOOST_CHECK_EQUAL(ts1.predicted(), ts2.predicted());
  BOOST_CHECK_EQUAL(ts1.predictedCovariance(), ts2.predictedCovariance());

  BOOST_CHECK_EQUAL(ts1.calibratedSourceLink(), ts2.calibratedSourceLink());
  BOOST_CHECK_EQUAL(ts1.calibrated(), ts2.calibrated());
  BOOST_CHECK_EQUAL(ts1.calibratedCovariance(), ts2.calibratedCovariance());
  BOOST_CHECK_EQUAL(ts1.calibratedSize(), ts2.calibratedSize());
  BOOST_CHECK_EQUAL(ts1.projector(), ts2.projector());

  BOOST_CHECK_EQUAL(ts1.jacobian(), ts2.jacobian());
  BOOST_CHECK_EQUAL(ts1.chi2(), ts2.chi2());              // always copied
  BOOST_CHECK_EQUAL(ts1.pathLength(), ts2.pathLength());  // always copied
  BOOST_CHECK_EQUAL(&ts1.referenceSurface(),
                    &ts2.referenceSurface());  // always copied
}

}  // namespace Test

}  // namespace Acts
