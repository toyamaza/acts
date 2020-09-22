// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "ActsExamples/Plugins/BField/ScalableBField.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include <random>
#include <stdexcept>

namespace {
template <typename TrackFinder>
struct TrackFinderFunctionImpl {
  TrackFinder trackFinder;

  TrackFinderFunctionImpl(TrackFinder&& f) : trackFinder(std::move(f)) {}

  ActsExamples::TrackFindingAlgorithm::TrackFinderResult operator()(
      const ActsExamples::SimSourceLinkContainer& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector>&
          options) const {
    return trackFinder.findTracks(sourceLinks, initialParameters, options);
  };
};
}  // namespace

ActsExamples::TrackFindingAlgorithm::TrackFinderFunction
ActsExamples::TrackFindingAlgorithm::makeTrackFinderFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    Options::BFieldVariant magneticField) {
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

  // unpack the magnetic field variant and instantiate the corresponding track
  // finder.
  return std::visit(
      [trackingGeometry](auto&& inputField) -> TrackFinderFunction {
        // each entry in the variant is already a shared_ptr
        // need ::element_type to get the real magnetic field type
        using InputMagneticField =
            typename std::decay_t<decltype(inputField)>::element_type;
        using MagneticField = Acts::SharedBField<InputMagneticField>;
        using Stepper = Acts::EigenStepper<MagneticField>;
        using Navigator = Acts::Navigator;
        using Propagator = Acts::Propagator<Stepper, Navigator>;
        using SourceLinkSelector = Acts::CKFSourceLinkSelector;
        using CKF =
            Acts::CombinatorialKalmanFilter<Propagator, Updater, Smoother,
                                            SourceLinkSelector>;

        // construct all components for the track finder
        MagneticField field(std::move(inputField));
        Stepper stepper(std::move(field));
        Navigator navigator(trackingGeometry);
        navigator.resolvePassive = false;
        navigator.resolveMaterial = true;
        navigator.resolveSensitive = true;
        Propagator propagator(std::move(stepper), std::move(navigator));
        CKF trackFinder(std::move(propagator));

        // build the track finder functions. owns the track finder object.
        return TrackFinderFunctionImpl<CKF>(std::move(trackFinder));
      },
      std::move(magneticField));
}
