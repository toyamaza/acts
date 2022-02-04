// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

template <typename E, typename R, typename A>
auto MultiEigenStepperLoop<E, R, A>::boundState(State& state,
                                                const Surface& surface,
                                                bool transportCov) const
    -> Result<BoundState> {
  if (numberComponents(state) == 1) {
    return SingleStepper::boundState(state.components.front().state, surface,
                                     transportCov);
  } else {  // Do the combination
    SmallVector<std::pair<double, BoundTrackParameters>> states;
    double accumulatedPathLength = 0.0;
    int failedBoundTransforms = 0;

    for (auto i = 0ul; i < numberComponents(state); ++i) {
      auto bs = SingleStepper::boundState(state.components[i].state, surface,
                                          transportCov);

      if (bs.ok()) {
        states.push_back(
            {state.components[i].weight, std::get<BoundTrackParameters>(*bs)});
        accumulatedPathLength +=
            std::get<double>(*bs) * state.components[i].weight;
      } else {
        failedBoundTransforms++;
      }
    }

    if (failedBoundTransforms > 0) {
      ACTS_ERROR("Multi component bound state: "
                 << failedBoundTransforms << " of " << numberComponents(state)
                 << " transforms failed");
    }

    if (states.size() == 0) {
      return MultiStepperError::AllComponentsConversionToBoundFailed;
    }

    // TODO At ATLAS, the final parameters seem to be computed with the mode of
    // the mixture. At the moment, we use the mean of the mixture here, but
    // there should be done a comparison sometimes in the future. This could
    // also be configurable maybe...
    const auto [params, cov] = detail::combineBoundGaussianMixture(
        states.begin(), states.end(), [&](const auto& wbs) {
          return std::tie(wbs.first, wbs.second.parameters(),
                          wbs.second.covariance());
        });

    return BoundState{BoundTrackParameters(surface.getSharedPtr(), params, cov),
                      Jacobian::Zero(), accumulatedPathLength};
  }
}

template <typename E, typename R, typename A>
auto MultiEigenStepperLoop<E, R, A>::curvilinearState(State& state,
                                                      bool transportCov) const
    -> CurvilinearState {
  if (numberComponents(state) == 1) {
    return SingleStepper::curvilinearState(state.components.front().state,
                                           transportCov);
  } else {
    Vector4 pos4 = Vector4::Zero();
    Vector3 dir = Vector3::Zero();
    ActsScalar qop = 0.0;
    BoundSymMatrix cov = BoundSymMatrix::Zero();
    ActsScalar pathLenth = 0.0;

    // TODO At ATLAS, the final parameters seem to be computed with the mode of
    // the mixture. At the moment, we use the mean of the mixture here, but
    // there should be done a comparison sometimes in the future. This could
    // also be configurable maybe...
    for (auto i = 0ul; i < numberComponents(state); ++i) {
      const auto [cp, jac, pl] = SingleStepper::curvilinearState(
          state.components[i].state, transportCov);

      pos4 += state.components[i].weight * cp.fourPosition(state.geoContext);
      dir += state.components[i].weight * cp.unitDirection();
      qop += state.components[i].weight * (cp.charge() / cp.absoluteMomentum());
      if (cp.covariance()) {
        cov += state.components[i].weight * *cp.covariance();
      }
      pathLenth += state.components[i].weight * pathLenth;
    }

    return CurvilinearState{CurvilinearTrackParameters(pos4, dir, qop, cov),
                            Jacobian::Zero(), pathLenth};
  }
}

template <typename E, typename R, typename A>
template <typename propagator_state_t>
Result<double> MultiEigenStepperLoop<E, R, A>::step(
    propagator_state_t& state) const {
  State& stepping = state.stepping;

  // Update step count
  stepping.steps++;

  // Check if we abort because of m_stepLimitAfterFirstComponentOnSurface
  if (stepping.stepCounterAfterFirstComponentOnSurface) {
    (*stepping.stepCounterAfterFirstComponentOnSurface)++;

    // If the limit is reached, remove all components which are not on a
    // surface, reweight the components, perform no step and return 0
    if (*stepping.stepCounterAfterFirstComponentOnSurface >=
        m_stepLimitAfterFirstComponentOnSurface) {
      auto& cmps = stepping.components;

      // It is not possible to remove components from the vector, since the
      // GSF actor relies on the fact that the ordering and number of
      // components does not change
      for (auto& cmp : cmps) {
        if (cmp.status != Intersection3D::Status::onSurface) {
          cmp.status = Intersection3D::Status::missed;
          cmp.weight = 0.0;
          cmp.state.pars.template segment<3>(eFreeDir0) = Vector3::Zero();
        }
      }

      // Reweight
      const auto sum_of_weights = std::accumulate(
          cmps.begin(), cmps.end(), ActsScalar{0},
          [](auto sum, const auto& cmp) { return sum + cmp.weight; });
      for (auto& cmp : cmps) {
        cmp.weight /= sum_of_weights;
      }

      ACTS_VERBOSE(
          "hit m_stepLimitAfterFirstComponentOnSurface, "
          "perform no step");

      stepping.stepCounterAfterFirstComponentOnSurface.reset();

      return 0.0;
    }
  }

  // Loop over all components and collect results in vector, write some
  // summary information to a stringstream
  SmallVector<Result<double>> results;
  std::stringstream ss;
  double accumulatedPathLength = 0.0;

  for (auto& component : stepping.components) {
    // We must also propagate missed components for the case that all
    // components miss the target and we need to re-target
    if (component.status == Intersection3D::Status::onSurface) {
      ss << "cmp skipped\t";
      continue;
    }

    using ThisSinglePropState =
        SinglePropState<SingleState, decltype(state.navigation),
                        decltype(state.options), decltype(state.geoContext)>;

    ThisSinglePropState single_state(component.state, state.navigation,
                                     state.options, state.geoContext);

    results.push_back(SingleStepper::step(single_state));

    if (results.back().ok()) {
      accumulatedPathLength += component.weight * *results.back();
      ss << *results.back() << "\t";
    } else {
      ss << "step error: " << results.back().error() << "\t";
    }
  }

  // Return no component was updated
  if (results.empty()) {
    return 0.0;
  }

  // Collect pointers to results which are ok, since Result is not copyable
  SmallVector<Result<double>*> ok_results;
  for (auto& res : results) {
    if (res.ok()) {
      ok_results.push_back(&res);
    }
  }

  // Return error if there is no ok result
  if (ok_results.empty()) {
    return MultiStepperError::AllComponentsSteppingError;
  }

  // Print the summary
  if (ok_results.size() == results.size()) {
    ACTS_VERBOSE("Performed steps: " << ss.str());
  } else {
    ACTS_WARNING("Performed steps with errors: " << ss.str());
  }

  // Return the weighted accumulated path length of all successful steps
  stepping.pathAccumulated += accumulatedPathLength;
  return accumulatedPathLength;
}
}  // namespace Acts
