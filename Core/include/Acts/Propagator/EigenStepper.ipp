// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"

template <typename B, typename E, typename A>
Acts::EigenStepper<B, E, A>::EigenStepper(B bField)
    : m_bField(std::move(bField)) {}

template <typename B, typename E, typename A>
void Acts::EigenStepper<B, E, A>::resetState(State& state,
                                             const BoundVector& boundParams,
                                             const BoundSymMatrix& cov,
                                             const Surface& surface,
                                             const NavigationDirection navDir,
                                             const double stepSize) const {
  // Update the stepping state
  update(state,
         detail::transformBoundToFreeParameters(surface, state.geoContext,
                                                boundParams),
         cov);
  state.navDir = navDir;
  state.stepSize = ConstrainedStep(stepSize);
  state.pathAccumulated = 0.;

  // Reinitialize the stepping jacobian
  surface.initJacobianToGlobal(state.geoContext, state.jacToGlobal,
                               position(state), direction(state), boundParams);
  state.jacobian = BoundMatrix::Identity();
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

template <typename B, typename E, typename A>
auto Acts::EigenStepper<B, E, A>::boundState(State& state,
                                             const Surface& surface) const
    -> BoundState {
  FreeVector parameters;
  parameters << state.pos[0], state.pos[1], state.pos[2], state.t, state.dir[0],
      state.dir[1], state.dir[2], state.q / state.p;
  return detail::boundState(state.geoContext, state.cov, state.jacobian,
                            state.jacTransport, state.derivative,
                            state.jacToGlobal, parameters, state.covTransport,
                            state.pathAccumulated, surface);
}

template <typename B, typename E, typename A>
auto Acts::EigenStepper<B, E, A>::curvilinearState(State& state) const
    -> CurvilinearState {
  FreeVector parameters;
  parameters << state.pos[0], state.pos[1], state.pos[2], state.t, state.dir[0],
      state.dir[1], state.dir[2], state.q / state.p;
  return detail::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, parameters, state.covTransport, state.pathAccumulated);
}

template <typename B, typename E, typename A>
void Acts::EigenStepper<B, E, A>::update(State& state,
                                         const FreeVector& parameters,
                                         const Covariance& covariance) const {
  state.pos = parameters.template segment<3>(eFreePos0);
  state.dir = parameters.template segment<3>(eFreeDir0).normalized();
  state.p = std::abs(1. / parameters[eFreeQOverP]);
  state.t = parameters[eFreeTime];

  state.cov = covariance;
}

template <typename B, typename E, typename A>
void Acts::EigenStepper<B, E, A>::update(State& state,
                                         const Vector3D& uposition,
                                         const Vector3D& udirection, double up,
                                         double time) const {
  state.pos = uposition;
  state.dir = udirection;
  state.p = up;
  state.t = time;
}

template <typename B, typename E, typename A>
void Acts::EigenStepper<B, E, A>::covarianceTransport(State& state) const {
  detail::covarianceTransport(state.cov, state.jacobian, state.jacTransport,
                              state.derivative, state.jacToGlobal, state.dir);
}

template <typename B, typename E, typename A>
void Acts::EigenStepper<B, E, A>::covarianceTransport(
    State& state, const Surface& surface) const {
  FreeVector parameters;
  parameters << state.pos[0], state.pos[1], state.pos[2], state.t, state.dir[0],
      state.dir[1], state.dir[2], state.q / state.p;
  detail::covarianceTransport(state.geoContext, state.cov, state.jacobian,
                              state.jacTransport, state.derivative,
                              state.jacToGlobal, parameters, surface);
}

template <typename B, typename E, typename A>
template <typename propagator_state_t>
Acts::Result<double> Acts::EigenStepper<B, E, A>::step(
    propagator_state_t& state) const {
  using namespace UnitLiterals;

  // Runge-Kutta integrator state
  auto& sd = state.stepping.stepData;
  double error_estimate = 0.;
  double h2, half_h;

  // First Runge-Kutta point (at current position)
  sd.B_first = getField(state.stepping, state.stepping.pos);
  if (!state.stepping.extension.validExtensionForStep(state, *this) ||
      !state.stepping.extension.k1(state, *this, sd.k1, sd.B_first, sd.kQoP)) {
    return 0.;
  }

  // The following functor starts to perform a Runge-Kutta step of a certain
  // size, going up to the point where it can return an estimate of the local
  // integration error. The results are stated in the local variables above,
  // allowing integration to continue once the error is deemed satisfactory
  const auto tryRungeKuttaStep = [&](const ConstrainedStep& h) -> bool {
    // State the square and half of the step size
    h2 = h * h;
    half_h = h * 0.5;

    // Second Runge-Kutta point
    const Vector3D pos1 =
        state.stepping.pos + half_h * state.stepping.dir + h2 * 0.125 * sd.k1;
    sd.B_middle = getField(state.stepping, pos1);
    if (!state.stepping.extension.k2(state, *this, sd.k2, sd.B_middle, sd.kQoP,
                                     half_h, sd.k1)) {
      return false;
    }

    // Third Runge-Kutta point
    if (!state.stepping.extension.k3(state, *this, sd.k3, sd.B_middle, sd.kQoP,
                                     half_h, sd.k2)) {
      return false;
    }

    // Last Runge-Kutta point
    const Vector3D pos2 =
        state.stepping.pos + h * state.stepping.dir + h2 * 0.5 * sd.k3;
    sd.B_last = getField(state.stepping, pos2);
    if (!state.stepping.extension.k4(state, *this, sd.k4, sd.B_last, sd.kQoP, h,
                                     sd.k3)) {
      return false;
    }

    // Compute and check the local integration error estimate
    error_estimate = std::max(
        h2 * ((sd.k1 - sd.k2 - sd.k3 + sd.k4).template lpNorm<1>() +
              std::abs(sd.kQoP[0] - sd.kQoP[1] - sd.kQoP[2] + sd.kQoP[3])),
        1e-20);
    return (error_estimate <= state.options.tolerance);
  };

  double stepSizeScaling = 1.;
  size_t nStepTrials = 0;
  // Select and adjust the appropriate Runge-Kutta step size as given
  // ATL-SOFT-PUB-2009-001
  while (!tryRungeKuttaStep(state.stepping.stepSize)) {
    stepSizeScaling =
        std::min(std::max(0.25, std::pow((state.options.tolerance /
                                          std::abs(2. * error_estimate)),
                                         0.25)),
                 4.);

    state.stepping.stepSize = state.stepping.stepSize * stepSizeScaling;

    // If step size becomes too small the particle remains at the initial
    // place
    if (state.stepping.stepSize * state.stepping.stepSize <
        state.options.stepSizeCutOff * state.options.stepSizeCutOff) {
      // Not moving due to too low momentum needs an aborter
      return EigenStepperError::StepSizeStalled;
    }

    // If the parameter is off track too much or given stepSize is not
    // appropriate
    if (nStepTrials > state.options.maxRungeKuttaStepTrials) {
      // Too many trials, have to abort
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
    nStepTrials++;
  }

  // use the adjusted step size
  const double h = state.stepping.stepSize;

  // When doing error propagation, update the associated Jacobian matrix
  if (state.stepping.covTransport) {
    // The step transport matrix in global coordinates
    FreeMatrix D;
    if (!state.stepping.extension.finalize(state, *this, h, D)) {
      return EigenStepperError::StepInvalid;
    }

    // for moment, only update the transport part
    state.stepping.jacTransport = D * state.stepping.jacTransport;
  } else {
    if (!state.stepping.extension.finalize(state, *this, h)) {
      return EigenStepperError::StepInvalid;
    }
  }

  // Update the track parameters according to the equations of motion
  state.stepping.pos +=
      h * state.stepping.dir + h2 / 6. * (sd.k1 + sd.k2 + sd.k3);
  state.stepping.dir += h / 6. * (sd.k1 + 2. * (sd.k2 + sd.k3) + sd.k4);
  state.stepping.dir /= state.stepping.dir.norm();
  if (state.stepping.covTransport) {
    state.stepping.derivative.template head<3>() = state.stepping.dir;
    state.stepping.derivative.template segment<3>(4) = sd.k4;
  }
  state.stepping.pathAccumulated += h;
  return h;
}
