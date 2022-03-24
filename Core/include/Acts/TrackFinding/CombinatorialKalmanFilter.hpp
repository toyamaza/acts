// This file is part of the Acts project.
//
// Copyright (C) 2016-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/TrackFinding/SourceLinkAccessorConcept.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/VoidKalmanComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <memory>
#include <unordered_map>

namespace Acts {

/// Track quality summary for one trajectory.
///
/// This could be used to decide if a track is to be recorded when the
/// filtering is done or to be terminated due to its bad quality
/// @todo: add other useful info, e.g. chi2
struct CombinatorialKalmanFilterTipState {
  // Number of passed sensitive surfaces
  size_t nSensitiveSurfaces = 0;
  // Number of track states
  size_t nStates = 0;
  // Number of (non-outlier) measurements
  size_t nMeasurements = 0;
  // Number of outliers
  size_t nOutliers = 0;
  // Number of holes
  size_t nHoles = 0;
};

/// Extension struct which holds the delegates to customize the CKF behavior
struct CombinatorialKalmanFilterExtensions {
  using candidate_container_t = std::vector<MultiTrajectory::TrackStateProxy>;
  using MeasurementSelector =
      Delegate<Result<std::pair<candidate_container_t::iterator,
                                candidate_container_t::iterator>>(
          candidate_container_t& trackStates, bool&, LoggerWrapper)>;
  using BranchStopper =
      Delegate<bool(const CombinatorialKalmanFilterTipState&)>;

  /// The Calibrator is a dedicated calibration algorithm that allows
  /// to calibrate measurements using track information, this could be
  /// e.g. sagging for wires, module deformations, etc.
  KalmanFitterExtensions::Calibrator calibrator;

  /// The updater incorporates measurement information into the track parameters
  KalmanFitterExtensions::Updater updater;

  /// The smoother back-propagates measurement information along the track
  KalmanFitterExtensions::Smoother smoother;

  /// The measurement selector is called during the filtering by the Actor.
  MeasurementSelector measurementSelector;

  BranchStopper branchStopper;

  /// Default constructor which connects the default void components
  CombinatorialKalmanFilterExtensions() {
    calibrator.connect<&voidKalmanCalibrator>();
    updater.connect<&voidKalmanUpdater>();
    smoother.connect<&voidKalmanSmoother>();
    branchStopper.connect<voidBranchStopper>();
    measurementSelector.connect<voidMeasurementSelector>();
  }

 private:
  /// Default measurement selector which will return all measurements
  /// @param candidates Measurement track state candidates
  /// @param isOutlier Output variable indicating whether the returned state is an outlier (unused)
  /// @param logger A logger instance
  static Result<
      std::pair<std::vector<MultiTrajectory::TrackStateProxy>::iterator,
                std::vector<MultiTrajectory::TrackStateProxy>::iterator>>
  voidMeasurementSelector(
      std::vector<MultiTrajectory::TrackStateProxy>& candidates,
      bool& isOutlier, LoggerWrapper logger) {
    (void)isOutlier;
    (void)logger;
    return std::pair{candidates.begin(), candidates.end()};
  };

  /// Default branch stopper which will never stop
  /// @param tipState The tip state to decide whether to stop (unused)
  /// @return false
  static bool voidBranchStopper(
      const CombinatorialKalmanFilterTipState& tipState) {
    (void)tipState;
    return false;
  }
};

/// Combined options for the combinatorial Kalman filter.
///
/// @tparam source_link_accessor_t Source link accessor type, should be
/// semiregular.
template <typename source_link_accessor_t>
struct CombinatorialKalmanFilterOptions {
  using SourceLinkAccessor = source_link_accessor_t;

  /// PropagatorOptions with context
  ///
  /// @param gctx The geometry context for this track finding/fitting
  /// @param mctx The magnetic context for this track finding/fitting
  /// @param cctx The calibration context for this track finding/fitting
  /// @param accessor_ The source link accessor
  /// @param extensions_ The extension struct
  /// @param logger_ The logger wrapper
  /// @param pOptions The plain propagator options
  /// @param rSurface The reference surface for the eventual track fitting to be
  /// expressed at
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param rSmoothing Whether to run smoothing to get fitted parameter
  CombinatorialKalmanFilterOptions(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      std::reference_wrapper<const CalibrationContext> cctx,
      SourceLinkAccessor accessor_,
      CombinatorialKalmanFilterExtensions extensions_, LoggerWrapper logger_,
      const PropagatorPlainOptions& pOptions, const Surface* rSurface = nullptr,
      bool mScattering = true, bool eLoss = true, bool rSmoothing = true)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        sourcelinkAccessor(std::move(accessor_)),
        extensions(std::move(extensions_)),
        propagatorPlainOptions(pOptions),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        smoothing(rSmoothing),
        logger(logger_) {}
  /// Contexts are required and the options must not be default-constructible.
  CombinatorialKalmanFilterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// The source link accessor
  SourceLinkAccessor sourcelinkAccessor;

  /// The filter extensions
  CombinatorialKalmanFilterExtensions extensions;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering.
  bool multipleScattering = true;

  /// Whether to consider energy loss.
  bool energyLoss = true;

  /// Whether to run smoothing to get fitted parameter
  bool smoothing = true;

  /// Logger instance
  LoggerWrapper logger;
};

struct CombinatorialKalmanFilterResult {
  // Fitted states that the actor has handled.
  MultiTrajectory fittedStates;

  // These is used internally to store candidate trackstates
  MultiTrajectory stateBuffer;
  std::vector<MultiTrajectory::TrackStateProxy> trackStateCandidates;

  // This is the indices of the 'tip' of the tracks stored in multitrajectory.
  // This correspond to the last measurment state in the multitrajectory.
  std::vector<size_t> lastMeasurementIndices;

  // This is the indices of the 'tip' of the tracks stored in multitrajectory.
  // This correspond to the last state in the multitrajectory.
  std::vector<size_t> lastTrackIndices;

  // The Parameters at the provided surface for separate tracks
  std::unordered_map<size_t, BoundTrackParameters> fittedParameters;

  // The indices of the 'tip' of the unfinished tracks
  std::vector<std::pair<size_t, CombinatorialKalmanFilterTipState>> activeTips;

  // The indices of track states and corresponding source links on different
  // surfaces
  std::unordered_map<const Surface*, std::unordered_map<size_t, size_t>>
      sourcelinkTips;

  // Indicator if filtering has been done
  bool filtered = false;

  // Indicator if smoothing has been done.
  bool smoothed = false;

  // The index for the current smoothing track
  size_t iSmoothed = 0;

  // Indicator if track finding has been done
  bool finished = false;

  Result<void> result{Result<void>::success()};
};

/// Combinatorial Kalman filter to find tracks.
///
///
/// @tparam propagator_t Type of the propagator
///
/// The CombinatorialKalmanFilter contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator
/// in order to initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update
/// and eventually the smoothing.  Updater, Smoother and Calibrator are
/// given to the Actor for further use:
/// - The Updater is the implemented kalman updater formalism, it
///   runs via a visitor pattern through the measurements.
/// - The Smoother is called at the end of the filtering (track finding) by the
/// Actor.
///
/// Measurements are not required to be ordered for the
/// CombinatorialKalmanFilter, measurement ordering needs to be figured out by
/// the navigation of the propagator.
///
/// The void components are provided mainly for unit testing.
template <typename propagator_t>
class CombinatorialKalmanFilter {
 public:
  /// Default constructor is deleted
  CombinatorialKalmanFilter() = delete;
  /// Constructor from arguments
  CombinatorialKalmanFilter(propagator_t pPropagator)
      : m_propagator(std::move(pPropagator)) {}

 private:
  using KalmanNavigator = typename propagator_t::Navigator;

  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// @brief Propagator Actor plugin for the CombinatorialKalmanFilter
  ///
  /// @tparam source_link_accessor_t The type of source link accessor
  /// @tparam parameters_t The type of parameters used for "local" paremeters.
  ///
  /// The CombinatorialKalmanFilter Actor does not rely on the measurements to
  /// be sorted along the track.
  template <typename source_link_accessor_t, typename parameters_t>
  class Actor {
   public:
    using TipState = CombinatorialKalmanFilterTipState;
    using BoundState = std::tuple<parameters_t, BoundMatrix, double>;
    using CurvilinearState =
        std::tuple<CurvilinearTrackParameters, BoundMatrix, double>;
    // The source link container type
    using SourceLinkContainer = typename source_link_accessor_t::Container;
    /// Broadcast the result_type
    using result_type = CombinatorialKalmanFilterResult;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Whether to run smoothing to get fitted parameter
    bool smoothing = true;

    /// @brief CombinatorialKalmanFilter actor operation
    ///
    /// @tparam propagator_state_t Type of the Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper is the stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void operator()(propagator_state_t& state, const stepper_t& stepper,
                    result_type& result) const {
      const auto& logger = state.options.logger;

      if (result.finished) {
        return;
      }

      ACTS_VERBOSE("CombinatorialKalmanFilter step");

      // Update:
      // - Waiting for a current surface
      auto surface = state.navigation.currentSurface;
      if (surface != nullptr and not result.filtered) {
        // There are three scenarios:
        // 1) The surface is in the measurement map
        // -> Select source links
        // -> Perform the kalman update for selected non-outlier source links
        // -> Add track states in multitrajectory. Multiple states mean branch
        // splitting.
        // -> Call branch stopper to justify each branch
        // -> If there is non-outlier state, update stepper information
        // 2) The surface is not in the measurement map but with material or is
        // an active surface
        // -> Add a hole or passive material state in multitrajectory
        // -> Call branch stopper to justify the branch
        // 3) The surface is neither in the measurement map nor with material
        // -> Do nothing
        ACTS_VERBOSE("Perform filter step");
        auto res = filter(surface, state, stepper, result);
        if (!res.ok()) {
          ACTS_ERROR("Error in filter: " << res.error());
          result.result = res.error();
        }
      }

      // Reset propagation state:
      // - When navigation breaks and there is stil active tip present after
      // recording&removing track tips on current surface
      if (state.navigation.navigationBreak and not result.filtered) {
        // Record the tips on current surface as trajectory entry indices
        // (taking advantage of fact that those tips are consecutive in list of
        // active tips) and remove those tips from active tips
        if (not result.activeTips.empty()) {
          // The last active tip
          const auto& lastActiveTip = result.activeTips.back().first;
          // Get the index of previous state
          const auto& iprevious =
              result.fittedStates.getTrackState(lastActiveTip).previous();
          // Find the track states which have the same previous state and remove
          // them from active tips
          while (not result.activeTips.empty()) {
            const auto& [currentTip, tipState] = result.activeTips.back();
            if (result.fittedStates.getTrackState(currentTip).previous() !=
                iprevious) {
              break;
            }
            // Record the tips if there are measurements on the track
            if (tipState.nMeasurements > 0) {
              ACTS_VERBOSE("Find track with entry index = "
                           << currentTip << " and there are nMeasurements = "
                           << tipState.nMeasurements
                           << ", nOutliers = " << tipState.nOutliers
                           << ", nHoles = " << tipState.nHoles << " on track");
              result.lastTrackIndices.emplace_back(currentTip);
              // Set the lastMeasurementIndex to the last measurement
              // to ignore the states after it in the rest of the algorithm
              auto lastMeasurementIndex = currentTip;
              auto lastMeasurementState =
                  result.fittedStates.getTrackState(lastMeasurementIndex);
              bool isMeasurement = lastMeasurementState.typeFlags().test(
                  TrackStateFlag::MeasurementFlag);
              while (!isMeasurement) {
                lastMeasurementIndex = lastMeasurementState.previous();
                lastMeasurementState =
                    result.fittedStates.getTrackState(lastMeasurementIndex);
                isMeasurement = lastMeasurementState.typeFlags().test(
                    TrackStateFlag::MeasurementFlag);
              }
              result.lastMeasurementIndices.emplace_back(lastMeasurementIndex);
            }
            // Remove the tip from list of active tips
            result.activeTips.erase(result.activeTips.end() - 1);
          }
        }
        // If no more active tip, done with filtering; Otherwise, reset
        // propagation state to track state at last tip of active tips
        if (result.activeTips.empty()) {
          ACTS_VERBOSE("Kalman filtering finds "
                       << result.lastTrackIndices.size() << " tracks");
          result.filtered = true;
        } else {
          ACTS_VERBOSE("Propagation jumps to branch with tip = "
                       << result.activeTips.back().first);
          reset(state, stepper, result);
        }
      }

      // Post-processing after filtering phase
      if (result.filtered) {
        // Return error if filtering finds no tracks
        if (result.lastTrackIndices.empty()) {
          result.result =
              Result<void>(CombinatorialKalmanFilterError::NoTrackFound);
        } else {
          if (not smoothing) {
            ACTS_VERBOSE("Finish Kalman filtering");
            // Remember that track finding is done
            result.finished = true;
          } else {
            // Iterate over the found tracks for smoothing and getting the
            // fitted parameter. This needs to be accomplished in different
            // propagation steps:
            // -> first run smoothing for found track indexed with iSmoothed
            if (not result.smoothed) {
              ACTS_VERBOSE(
                  "Finalize/run smoothing for track with last measurement "
                  "index = "
                  << result.lastMeasurementIndices.at(result.iSmoothed));
              // --> Search the starting state to run the smoothing
              // --> Call the smoothing
              // --> Set a stop condition when all track states have been
              // handled
              auto res = finalize(state, stepper, result);
              if (!res.ok()) {
                ACTS_ERROR("Error in finalize: " << res.error());
                result.result = res.error();
              }
              result.smoothed = true;
            }
            // -> then progress to target/reference surface and built the final
            // track parameters for found track indexed with iSmoothed
            if (result.smoothed and
                targetReached(state, stepper, *targetSurface)) {
              ACTS_VERBOSE(
                  "Completing the track with last measurement index = "
                  << result.lastMeasurementIndices.at(result.iSmoothed));
              // Transport & bind the parameter to the final surface
              auto res = stepper.boundState(state.stepping, *targetSurface);
              if (!res.ok()) {
                ACTS_ERROR("Error in finalize: " << res.error());
                result.result = res.error();
                return;
              }

              auto fittedState = *res;
              // Assign the fitted parameters
              result.fittedParameters.emplace(
                  result.lastMeasurementIndices.at(result.iSmoothed),
                  std::get<BoundTrackParameters>(fittedState));
              // If there are more trajectories to handle:
              // -> set the targetReached status to false
              // -> set the smoothed status to false
              // -> update the index of track to be smoothed
              if (result.iSmoothed < result.lastMeasurementIndices.size() - 1) {
                state.navigation.targetReached = false;
                result.smoothed = false;
                result.iSmoothed++;
                // Reverse navigation direction to start targeting for the rest
                // tracks
                state.stepping.navDir =
                    (state.stepping.navDir == backward) ? forward : backward;
                // To avoid meaningless navigation target call
                state.stepping.stepSize =
                    ConstrainedStep(state.stepping.navDir *
                                    std::abs(state.options.maxStepSize));
              } else {
                ACTS_VERBOSE("Finish Kalman filtering and smoothing");
                // Remember that track finding is done
                result.finished = true;
              }
            }
          }  // if run smoothing
        }    // if there are found tracks
      }      // if filtering is done
    }

    /// @brief Kalman actor operation : reset propagation
    ///
    /// @tparam propagator_state_t Type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper is the stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void reset(propagator_state_t& state, stepper_t& stepper,
               result_type& result) const {
      auto currentState =
          result.fittedStates.getTrackState(result.activeTips.back().first);

      // Update the stepping state
      stepper.resetState(state.stepping, currentState.filtered(),
                         currentState.filteredCovariance(),
                         currentState.referenceSurface(), state.stepping.navDir,
                         state.options.maxStepSize);

      // Reset the navigation state
      // Set targetSurface to nullptr for forward filtering; it's only needed
      // after smoothing
      state.navigation.reset(state.geoContext, stepper.position(state.stepping),
                             stepper.direction(state.stepping),
                             state.stepping.navDir,
                             &currentState.referenceSurface(), nullptr);

      // No Kalman filtering for the starting surface, but still need
      // to consider the material effects here
      materialInteractor(state.navigation.currentSurface, state, stepper);
    }

    /// @brief CombinatorialKalmanFilter actor operation :
    /// - filtering for all measurement(s) on surface
    /// - store selected track states in multiTrajectory
    /// - update propagator state to the (last) selected track state
    ///
    /// @tparam propagator_state_t Type of the Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    Result<void> filter(const Surface* surface, propagator_state_t& state,
                        const stepper_t& stepper, result_type& result) const {
      const auto& logger = state.options.logger;
      // Initialize the number of branches on current surface
      size_t nBranchesOnSurface = 0;

      // Count the number of source links on the surface
      size_t nSourcelinks = m_sourcelinkAccessor.count(surface->geometryId());
      if (nSourcelinks > 0) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface " << surface->geometryId()
                                            << " detected.");

        // Transport the covariance to the surface
        stepper.transportCovarianceToBound(state.stepping, *surface);

        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, preUpdate);

        // Bind the transported state to the current surface
        auto boundStateRes =
            stepper.boundState(state.stepping, *surface, false);
        if (!boundStateRes.ok()) {
          return boundStateRes.error();
        }
        auto boundState = *boundStateRes;

        // Retrieve the previous tip and its state
        // The states created on this surface will have the common previous tip
        size_t prevTip = SIZE_MAX;
        TipState prevTipState;
        if (not result.activeTips.empty()) {
          prevTip = result.activeTips.back().first;
          prevTipState = result.activeTips.back().second;
          // New state is to be added. Remove the last tip from active tips
          result.activeTips.erase(result.activeTips.end() - 1);
        }

        // Create trackstates for all source links (will be filtered later)
        // Results are stored in result => no return value
        createSourceLinkTrackStates(state.geoContext, surface, result,
                                    boundState, prevTip);

        // Invoke the measurement selector to select compatible measurements
        // with the predicted track parameter.
        // It can modify the trackStateCandidates vector, and will return a pair
        // of iterators marking the range of accepted measurements (track
        // states)
        bool isOutlier = false;
        auto selectorResult = m_extensions.measurementSelector(
            result.trackStateCandidates, isOutlier, logger);

        if (!selectorResult.ok()) {
          ACTS_ERROR("Selection of calibrated measurements failed: "
                     << selectorResult.error());
          return selectorResult.error();
        }
        auto selectedTrackStateRange = *selectorResult;

        auto procRes = processSelectedTrackStates(
            state.geoContext, selectedTrackStateRange.first,
            selectedTrackStateRange.second, result, isOutlier, prevTipState,
            nBranchesOnSurface, logger);

        if (!procRes.ok()) {
          ACTS_ERROR(
              "Processing of selected track states failed: " << procRes.error())
          return procRes.error();
        }

        if (nBranchesOnSurface > 0 and not isOutlier) {
          // If there are measurement track states on this surface
          ACTS_VERBOSE("Filtering step successful with " << nBranchesOnSurface
                                                         << " branches");
          // Update stepping state using filtered parameters of last track
          // state on this surface
          auto ts =
              result.fittedStates.getTrackState(result.activeTips.back().first);
          stepper.update(state.stepping,
                         MultiTrajectoryHelpers::freeFiltered(
                             state.options.geoContext, ts),
                         ts.filtered(), ts.filteredCovariance(), *surface);
          ACTS_VERBOSE("Stepping state is updated with filtered parameter:");
          ACTS_VERBOSE("-> " << ts.filtered().transpose()
                             << " of track state with tip = "
                             << result.activeTips.back().first);
        }
        // Update state and stepper with post material effects
        materialInteractor(surface, state, stepper, postUpdate);
      } else if (surface->associatedDetectorElement() != nullptr ||
                 surface->surfaceMaterial() != nullptr) {
        // No splitting on the surface without source links. Set it to one
        // first, but could be changed later
        nBranchesOnSurface = 1;

        // Retrieve the previous tip and its state
        size_t prevTip = SIZE_MAX;
        TipState tipState;
        if (not result.activeTips.empty()) {
          prevTip = result.activeTips.back().first;
          tipState = result.activeTips.back().second;
        }

        // The surface could be either sensitive or passive
        bool isSensitive = (surface->associatedDetectorElement() != nullptr);
        bool isMaterial = (surface->surfaceMaterial() != nullptr);
        std::string type = isSensitive ? "sensitive" : "passive";
        ACTS_VERBOSE("Detected " << type
                                 << " surface: " << surface->geometryId());
        if (isSensitive) {
          // Increment of number of passed sensitive surfaces
          tipState.nSensitiveSurfaces++;
        }
        // Add state if there is already measurement detected on this branch
        // For in-sensitive surface, only add state when smoothing is
        // required
        bool createState = false;
        if (smoothing) {
          createState = (tipState.nMeasurements > 0 or isMaterial);
        } else {
          createState = (tipState.nMeasurements > 0 and isSensitive);
        }
        if (createState) {
          // New state is to be added. Remove the last tip from active tips now
          if (not result.activeTips.empty()) {
            result.activeTips.erase(result.activeTips.end() - 1);
          }
          // No source links on surface, add either hole or passive material
          // TrackState. No storage allocation for uncalibrated/calibrated
          // measurement and filtered parameter
          auto stateMask =
              ~(TrackStatePropMask::Uncalibrated |
                TrackStatePropMask::Calibrated | TrackStatePropMask::Filtered);

          // Increment of number of processed states
          tipState.nStates++;
          size_t currentTip = SIZE_MAX;
          if (isSensitive) {
            // Incremet of number of holes
            tipState.nHoles++;
          }

          // Transport & bind the state to the current surface
          auto res = stepper.boundState(state.stepping, *surface);
          if (!res.ok()) {
            ACTS_ERROR("Error in filter: " << res.error());
            return res.error();
          }
          const auto boundState = *res;
          // Add a hole or material track state to the multitrajectory
          currentTip = addNonSourcelinkState(stateMask, boundState, result,
                                             isSensitive, prevTip, logger);

          // Check the branch
          if (not m_extensions.branchStopper(tipState)) {
            // Remember the active tip and its state
            result.activeTips.emplace_back(std::move(currentTip),
                                           std::move(tipState));
          } else {
            // No branch on this surface
            nBranchesOnSurface = 0;
          }
        }
        if (surface->surfaceMaterial() != nullptr) {
          // Update state and stepper with material effects
          materialInteractor(surface, state, stepper, fullUpdate);
        }
      } else {
        // Neither measurement nor material on surface, this branch is still
        // valid. Count the branch on current surface
        nBranchesOnSurface = 1;
      }

      // Reset current tip if there is no branch on current surface
      if (nBranchesOnSurface == 0) {
        ACTS_DEBUG("Branch on surface " << surface->geometryId()
                                        << " is stopped");
        if (not result.activeTips.empty()) {
          ACTS_VERBOSE("Propagation jumps to branch with tip = "
                       << result.activeTips.back().first);
          reset(state, stepper, result);
        } else {
          ACTS_VERBOSE("Stop Kalman filtering with "
                       << result.lastMeasurementIndices.size()
                       << " found tracks");
          result.filtered = true;
        }
      }

      return Result<void>::success();
    }

    /// Create and fill track states for all source links
    /// @param gctx The current geometry context
    /// @param surface The surface currently being processed
    /// @param result Reference to the result struct of the actor
    /// @param boundState Bound state from the propagation on this surface
    /// @param prevTip Index pointing at previous trajectory state (i.e. tip)
    void createSourceLinkTrackStates(const Acts::GeometryContext& gctx,
                                     const Surface* surface,
                                     result_type& result,
                                     const BoundState& boundState,
                                     size_t prevTip) const {
      const auto& [boundParams, jacobian, pathLength] = boundState;

      // Get all source links on the surface
      auto [lower_it, upper_it] =
          m_sourcelinkAccessor.range(surface->geometryId());

      result.trackStateCandidates.clear();
      result.trackStateCandidates.reserve(std::distance(lower_it, upper_it));

      result.stateBuffer.clear();

      using PM = TrackStatePropMask;

      // Calibrate all the source links on the surface since the selection has
      // to be done based on calibrated measurement
      for (auto it = lower_it; it != upper_it; ++it) {
        // get the source link
        const auto& sourceLink = m_sourcelinkAccessor.at(it);

        // prepare the track state
        PM mask =
            PM::Predicted | PM::Jacobian | PM::Uncalibrated | PM::Calibrated;

        if (it != lower_it) {
          // not the first TrackState, only need uncalibrated and calibrated
          mask = PM::Uncalibrated | PM::Calibrated;
        }

        size_t tsi = result.stateBuffer.addTrackState(mask, prevTip);
        // CAREFUL! This trackstate has a previous index that is not in this
        // MultiTrajectory Visiting brackwards from this track state will
        // fail!
        auto ts = result.stateBuffer.getTrackState(tsi);

        if (it == lower_it) {
          // only set these for first
          ts.predicted() = boundParams.parameters();
          if (boundParams.covariance()) {
            ts.predictedCovariance() = *boundParams.covariance();
          }
          ts.jacobian() = jacobian;
        } else {
          // subsequent track states can reuse
          auto& first = result.trackStateCandidates.front();
          ts.data().ipredicted = first.data().ipredicted;
          ts.data().ijacobian = first.data().ijacobian;
        }

        ts.pathLength() = pathLength;

        ts.setReferenceSurface(boundParams.referenceSurface().getSharedPtr());

        ts.setUncalibrated(sourceLink);

        // now calibrate the track state
        m_extensions.calibrator(gctx, ts);

        result.trackStateCandidates.push_back(ts);
      }
    }

    /// Handle the list of selected track states
    /// @param gctx The current geometry context
    /// @param begin The start iterator for selected track states
    /// @param end The end iterator for selected track states
    /// @param result Reference to the actor result struct
    /// @param isOutlier If this track state is a single outlier one
    /// @param prevTipState Tip state prior to this surface
    /// @param [in,out] nBranchesOnSurface Number of branches on surface, will be updated
    /// @param logger A logging instance
    Result<void> processSelectedTrackStates(
        const Acts::GeometryContext& gctx,
        std::vector<MultiTrajectory::TrackStateProxy>::const_iterator begin,
        std::vector<MultiTrajectory::TrackStateProxy>::const_iterator end,
        result_type& result, bool isOutlier, const TipState& prevTipState,
        size_t& nBranchesOnSurface, LoggerWrapper logger) const {
      using PM = TrackStatePropMask;

      std::optional<MultiTrajectory::TrackStateProxy> firstTrackState{
          std::nullopt};
      for (auto it = begin; it != end; ++it) {
        auto& candidateTrackState = *it;

        PM mask = PM::All;

        if (it != begin) {
          // subsequent track states don't need storage for these
          mask = ~PM::Predicted & ~PM::Jacobian;
        }

        if (isOutlier) {
          mask &= ~PM::Filtered;  // outlier won't have separate filtered
                                  // parameters
        }

        // copy this trackstate into fitted states MultiTrajectory
        MultiTrajectory::TrackStateProxy trackState =
            result.fittedStates.getTrackState(result.fittedStates.addTrackState(
                mask, candidateTrackState.previous()));

        if (it != begin) {
          // assign indices pointing to first track state
          trackState.data().ipredicted = firstTrackState->data().ipredicted;
          trackState.data().ijacobian = firstTrackState->data().ijacobian;
        } else {
          firstTrackState = trackState;
        }

        // either copy ALL or everything except for predicted and jacobian
        trackState.copyFrom(candidateTrackState, mask, false);

        auto& typeFlags = trackState.typeFlags();
        if (trackState.referenceSurface().surfaceMaterial() != nullptr) {
          typeFlags.set(TrackStateFlag::MaterialFlag);
        }
        typeFlags.set(TrackStateFlag::ParameterFlag);

        // Inherit the tip state from the previous and will be updated
        // later
        TipState tipState = prevTipState;
        size_t currentTip = trackState.index();

        // Increment of number of processedState and passed sensitive surfaces
        tipState.nSensitiveSurfaces++;
        tipState.nStates++;

        if (isOutlier) {
          ACTS_VERBOSE(
              "Creating outlier track state with tip = " << currentTip);
          // Set the outlier flag
          typeFlags.set(TrackStateFlag::OutlierFlag);
          // Increment number of outliers
          tipState.nOutliers++;
          // No Kalman update for outlier
          // Set the filtered parameter index to be the same with predicted
          // parameter
          trackState.data().ifiltered = trackState.data().ipredicted;

        } else {
          // Kalman update
          auto updateRes =
              m_extensions.updater(gctx, trackState, forward, getDummyLogger());
          if (!updateRes.ok()) {
            ACTS_ERROR("Update step failed: " << updateRes.error());
            return updateRes.error();
          }
          ACTS_VERBOSE(
              "Creating measurement track state with tip = " << currentTip);
          // Set the measurement flag
          typeFlags.set(TrackStateFlag::MeasurementFlag);
          // Increment number of measurements
          tipState.nMeasurements++;
        }

        // Check if need to stop this branch
        if (not m_extensions.branchStopper(tipState)) {
          // Put tipstate back into active tips to continue with it
          result.activeTips.emplace_back(std::move(currentTip),
                                         std::move(tipState));
          // Record the number of branches on surface
          nBranchesOnSurface++;
        }
      }
      return Result<void>::success();
    }

    /// @brief CombinatorialKalmanFilter actor operation : add hole or material track state
    ///
    /// @param stateMask The bitmask that instructs which components to allocate
    /// @param boundState The bound state on current surface
    /// @param result is the mutable result state object
    /// and which to leave invalid
    /// @param isSensitive The surface is sensitive or passive
    /// @param prevTip The index of the previous state
    /// @param logger The logger wrapper
    ///
    /// @return The tip of added state
    size_t addNonSourcelinkState(
        const TrackStatePropMask& stateMask, const BoundState& boundState,
        result_type& result, bool isSensitive, size_t prevTip = SIZE_MAX,
        LoggerWrapper logger = getDummyLogger()) const {
      // Add a track state
      auto currentTip = result.fittedStates.addTrackState(stateMask, prevTip);
      if (isSensitive) {
        ACTS_VERBOSE("Creating Hole track state with tip = " << currentTip);
      } else {
        ACTS_VERBOSE("Creating Material track state with tip = " << currentTip);
      }
      // now get track state proxy back
      auto trackStateProxy = result.fittedStates.getTrackState(currentTip);

      const auto& [boundParams, jacobian, pathLength] = boundState;
      // Fill the track state
      trackStateProxy.predicted() = boundParams.parameters();
      if (boundParams.covariance().has_value()) {
        trackStateProxy.predictedCovariance() = *boundParams.covariance();
      }
      trackStateProxy.jacobian() = jacobian;
      trackStateProxy.pathLength() = pathLength;
      // Set the surface
      trackStateProxy.setReferenceSurface(
          boundParams.referenceSurface().getSharedPtr());
      // Set the filtered parameter index to be the same with predicted
      // parameter

      // Set the track state flags
      auto& typeFlags = trackStateProxy.typeFlags();
      if (trackStateProxy.referenceSurface().surfaceMaterial() != nullptr) {
        typeFlags.set(TrackStateFlag::MaterialFlag);
      }
      typeFlags.set(TrackStateFlag::ParameterFlag);
      if (isSensitive) {
        typeFlags.set(TrackStateFlag::HoleFlag);
      }

      trackStateProxy.data().ifiltered = trackStateProxy.data().ipredicted;

      return currentTip;
    }

    /// @brief CombinatorialKalmanFilter actor operation : material interaction
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the material interaction happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param updateStage The materal update stage
    ///
    template <typename propagator_state_t, typename stepper_t>
    void materialInteractor(
        const Surface* surface, propagator_state_t& state, stepper_t& stepper,
        const MaterialUpdateStage& updateStage = fullUpdate) const {
      const auto& logger = state.options.logger;
      // Indicator if having material
      bool hasMaterial = false;

      if (surface and surface->surfaceMaterial()) {
        // Prepare relevant input particle properties
        detail::PointwiseMaterialInteraction interaction(surface, state,
                                                         stepper);
        // Evaluate the material properties
        if (interaction.evaluateMaterialSlab(state, updateStage)) {
          // Surface has material at this stage
          hasMaterial = true;

          // Evaluate the material effects
          interaction.evaluatePointwiseMaterialInteraction(multipleScattering,
                                                           energyLoss);

          // Screen out material effects info
          ACTS_VERBOSE("Material effects on surface: "
                       << surface->geometryId()
                       << " at update stage: " << updateStage << " are :");
          ACTS_VERBOSE("eLoss = "
                       << interaction.Eloss << ", "
                       << "variancePhi = " << interaction.variancePhi << ", "
                       << "varianceTheta = " << interaction.varianceTheta
                       << ", "
                       << "varianceQoverP = " << interaction.varianceQoverP);

          // Update the state and stepper with material effects
          interaction.updateState(state, stepper, addNoise);
        }
      }

      if (not hasMaterial) {
        // Screen out message
        ACTS_VERBOSE("No material effects on surface: " << surface->geometryId()
                                                        << " at update stage: "
                                                        << updateStage);
      }
    }

    /// @brief Kalman actor operation : finalize
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    Result<void> finalize(propagator_state_t& state, const stepper_t& stepper,
                          result_type& result) const {
      const auto& logger = state.options.logger;
      // The measurement tip of the track being smoothed
      const auto& lastMeasurementIndex =
          result.lastMeasurementIndices.at(result.iSmoothed);

      // Get the indices of the first states (can be either a measurement or
      // material);
      size_t firstStateIndex = lastMeasurementIndex;
      // Count track states to be smoothed
      size_t nStates = 0;
      result.fittedStates.applyBackwards(lastMeasurementIndex, [&](auto st) {
        bool isMeasurement =
            st.typeFlags().test(TrackStateFlag::MeasurementFlag);
        bool isMaterial = st.typeFlags().test(TrackStateFlag::MaterialFlag);
        if (isMeasurement || isMaterial) {
          firstStateIndex = st.index();
        }
        nStates++;
      });
      // Return error if the track has no measurement states (but this should
      // not happen)
      if (nStates == 0) {
        ACTS_ERROR("Smoothing for a track without measurements.");
        return CombinatorialKalmanFilterError::SmoothFailed;
      }
      // Screen output for debugging
      ACTS_VERBOSE("Apply smoothing on " << nStates
                                         << " filtered track states.");
      // Smooth the track states
      auto smoothRes =
          m_extensions.smoother(state.geoContext, result.fittedStates,
                                lastMeasurementIndex, getDummyLogger());
      if (!smoothRes.ok()) {
        ACTS_ERROR("Smoothing step failed: " << smoothRes.error());
        return smoothRes.error();
      }

      // Return in case no target surface
      if (targetSurface == nullptr) {
        return Result<void>::success();
      }

      // Obtain the smoothed parameters at first/last measurement state.
      // The first state can also be a material state
      auto firstCreatedState =
          result.fittedStates.getTrackState(firstStateIndex);
      auto lastCreatedMeasurement =
          result.fittedStates.getTrackState(lastMeasurementIndex);

      // Lambda to get the intersection of the free params on the target surface
      auto target = [&](const FreeVector& freeVector) -> SurfaceIntersection {
        return targetSurface->intersect(
            state.geoContext, freeVector.segment<3>(eFreePos0),
            state.stepping.navDir * freeVector.segment<3>(eFreeDir0), true);
      };

      // The smoothed free params at the first/last measurement state
      auto firstParams = MultiTrajectoryHelpers::freeSmoothed(
          state.options.geoContext, firstCreatedState);
      auto lastParams = MultiTrajectoryHelpers::freeSmoothed(
          state.options.geoContext, lastCreatedMeasurement);
      // Get the intersections of the smoothed free parameters with the target
      // surface
      const auto firstIntersection = target(firstParams);
      const auto lastIntersection = target(lastParams);

      // Update the stepping parameters - in order to progress to destination.
      // At the same time, reverse navigation direction for further
      // stepping if necessary.
      // @note The stepping parameters is updated to the smoothed parameters at
      // either the first measurement state or the last measurement state. It
      // assumes the target surface is not within the first and the last
      // smoothed measurement state. Also, whether the intersection is on
      // surface is not checked here.
      bool reverseDirection = false;
      bool closerToFirstCreatedState =
          (std::abs(firstIntersection.intersection.pathLength) <=
           std::abs(lastIntersection.intersection.pathLength));
      if (closerToFirstCreatedState) {
        stepper.update(state.stepping, firstParams,
                       firstCreatedState.smoothed(),
                       firstCreatedState.smoothedCovariance(),
                       firstCreatedState.referenceSurface());
        reverseDirection = (firstIntersection.intersection.pathLength < 0);
      } else {
        stepper.update(state.stepping, lastParams,
                       lastCreatedMeasurement.smoothed(),
                       lastCreatedMeasurement.smoothedCovariance(),
                       lastCreatedMeasurement.referenceSurface());
        reverseDirection = (lastIntersection.intersection.pathLength < 0);
      }
      const auto& surface = closerToFirstCreatedState
                                ? firstCreatedState.referenceSurface()
                                : lastCreatedMeasurement.referenceSurface();
      ACTS_VERBOSE(
          "Smoothing successful, updating stepping state to smoothed "
          "parameters at surface "
          << surface.geometryId() << ". Prepared to reach the target surface.");

      // Reverse the navigation direction if necessary
      if (reverseDirection) {
        ACTS_VERBOSE(
            "Reverse navigation direction after smoothing for reaching the "
            "target surface");
        state.stepping.navDir =
            (state.stepping.navDir == forward) ? backward : forward;
      }
      // Reinitialize the stepping jacobian
      state.stepping.jacobian = BoundMatrix::Identity();
      state.stepping.jacTransport = FreeMatrix::Identity();
      state.stepping.derivative = FreeVector::Zero();
      // Reset the step size
      state.stepping.stepSize = ConstrainedStep(
          state.stepping.navDir * std::abs(state.options.maxStepSize));
      // Set accumulatd path to zero before targeting surface
      state.stepping.pathAccumulated = 0.;

      return Result<void>::success();
    }

    CombinatorialKalmanFilterExtensions m_extensions;

    /// The source link accesor
    source_link_accessor_t m_sourcelinkAccessor;

    /// The Surface being targeted
    SurfaceReached targetReached;
  };

  template <typename source_link_accessor_t, typename parameters_t>
  class Aborter {
   public:
    /// Broadcast the result_type
    using action_type = Actor<source_link_accessor_t, parameters_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename result_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const result_t& result) const {
      if (!result.result.ok() or result.finished) {
        return true;
      }
      return false;
    }
  };

 public:
  /// Combinatorial Kalman Filter implementation, calls the the Kalman filter
  /// and smoother
  ///
  /// @tparam source_link_accessor_t Type of the source link accessor
  /// @tparam start_parameters_container_t Type of the initial parameters
  /// container
  /// @tparam calibrator_t Type of the source link calibrator
  /// @tparam measurement_selector_t Type of the measurement selector
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param initialParameters The initial track parameters
  /// @param tfOptions CombinatorialKalmanFilterOptions steering the track
  /// finding
  /// @note The input measurements are given in the form of @c SourceLinks.
  /// It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the track finding.
  ///
  /// @return a container of track finding result for all the initial track
  /// parameters
  template <typename source_link_accessor_t,
            typename start_parameters_container_t,
            typename parameters_t = BoundTrackParameters>
  std::vector<Result<CombinatorialKalmanFilterResult>> findTracks(
      const typename source_link_accessor_t::Container& sourcelinks,
      const start_parameters_container_t& initialParameters,
      const CombinatorialKalmanFilterOptions<source_link_accessor_t>& tfOptions)
      const {
    static_assert(
        SourceLinkAccessorConcept<source_link_accessor_t>,
        "The source link accessor does not fullfill SourceLinkAccessorConcept");
    static_assert(
        std::is_same_v<GeometryIdentifier,
                       typename source_link_accessor_t::Key>,
        "The source link container does not have GeometryIdentifier as the key "
        "type");

    const auto& logger = tfOptions.logger;

    ACTS_VERBOSE("Preparing " << sourcelinks.size() << " input measurements");

    // Create the ActionList and AbortList
    using CombinatorialKalmanFilterAborter =
        Aborter<source_link_accessor_t, parameters_t>;
    using CombinatorialKalmanFilterActor =
        Actor<source_link_accessor_t, parameters_t>;
    using Actors = ActionList<CombinatorialKalmanFilterActor>;
    using Aborters = AbortList<CombinatorialKalmanFilterAborter>;

    // Create relevant options for the propagation options
    PropagatorOptions<Actors, Aborters> propOptions(
        tfOptions.geoContext, tfOptions.magFieldContext, tfOptions.logger);

    // Set the trivial propagator options
    propOptions.setPlainOptions(tfOptions.propagatorPlainOptions);

    // Catch the actor
    auto& combKalmanActor =
        propOptions.actionList.template get<CombinatorialKalmanFilterActor>();
    combKalmanActor.targetSurface = tfOptions.referenceSurface;
    combKalmanActor.multipleScattering = tfOptions.multipleScattering;
    combKalmanActor.energyLoss = tfOptions.energyLoss;
    combKalmanActor.smoothing = tfOptions.smoothing;

    // copy source link accessor, calibrator and measurement selector
    combKalmanActor.m_sourcelinkAccessor = tfOptions.sourcelinkAccessor;
    // set the pointer to the source links
    combKalmanActor.m_sourcelinkAccessor.container = &sourcelinks;
    combKalmanActor.m_extensions = tfOptions.extensions;

    // Run the CombinatorialKalmanFilter.
    // @todo The same target surface is used for all the initial track
    // parameters, which is not necessarily the case.
    std::vector<Result<CombinatorialKalmanFilterResult>> ckfResults;
    ckfResults.reserve(initialParameters.size());
    // Loop over all initial track parameters. Return the results for all
    // initial track parameters including those failed ones.
    for (size_t iseed = 0; iseed < initialParameters.size(); ++iseed) {
      const auto& sParameters = initialParameters[iseed];
      auto result = m_propagator.template propagate(sParameters, propOptions);

      if (!result.ok()) {
        ACTS_ERROR("Propapation failed: "
                   << result.error() << " " << result.error().message()
                   << " with the initial parameters " << iseed << " : \n"
                   << sParameters.parameters());
        // Emplace back the failed result
        ckfResults.emplace_back(result.error());
        continue;
      }

      const auto& propRes = *result;

      /// Get the result of the CombinatorialKalmanFilter
      auto combKalmanResult =
          propRes.template get<CombinatorialKalmanFilterResult>();

      /// The propagation could already reach max step size
      /// before the track finding is finished during two phases:
      // -> filtering for track finding;
      // -> surface targeting to get fitted parameters at target surface.
      // This is regarded as a failure.
      // @TODO: Implement distinguishment between the above two cases if
      // necessary
      if (combKalmanResult.result.ok() and not combKalmanResult.finished) {
        combKalmanResult.result = Result<void>(
            CombinatorialKalmanFilterError::PropagationReachesMaxSteps);
      }

      if (!combKalmanResult.result.ok()) {
        ACTS_ERROR("CombinatorialKalmanFilter failed: "
                   << combKalmanResult.result.error() << " "
                   << combKalmanResult.result.error().message()
                   << " with the initial parameters " << iseed << " : \n"
                   << sParameters.parameters());
        // Emplace back the failed result
        ckfResults.emplace_back(combKalmanResult.result.error());
        continue;
      }

      // Emplace back the successful result
      ckfResults.emplace_back(combKalmanResult);
    }

    return ckfResults;
  }

};  // namespace Acts

}  // namespace Acts
