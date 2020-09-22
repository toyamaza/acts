// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

ActsExamples::TrackFindingAlgorithm::TrackFindingAlgorithm(
    Config cfg, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("TrackFindingAlgorithm", level),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing input source links collection");
  }
  if (m_cfg.inputInitialTrackParameters.empty()) {
    throw std::invalid_argument(
        "Missing input initial track parameters collection");
  }
  if (m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument("Missing output trajectories collection");
  }
}

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  const auto sourceLinks =
      ctx.eventStore.get<SimSourceLinkContainer>(m_cfg.inputSourceLinks);
  const auto initialParameters = ctx.eventStore.get<TrackParametersContainer>(
      m_cfg.inputInitialTrackParameters);

  // Prepare the output data with MultiTrajectory
  TrajectoryContainer trajectories;
  trajectories.reserve(initialParameters.size());

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3D{0., 0., 0.});

  Acts::PropagatorPlainOptions pOptions;
  pOptions.maxSteps = 10000;

  // Set the CombinatorialKalmanFilter options
  ActsExamples::TrackFindingAlgorithm::CKFOptions ckfOptions(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
      m_cfg.sourcelinkSelectorCfg, Acts::LoggerWrapper{logger()}, pOptions,
      &(*pSurface));

  // Perform the track finding for each starting parameter
  // @TODO: use seeds from track seeding algorithm as starting parameter
  for (std::size_t iseed = 0; iseed < initialParameters.size(); ++iseed) {
    const auto& initialParams = initialParameters[iseed];

    ACTS_DEBUG("Invoke track finding seeded by truth particle " << iseed);
    auto result = m_cfg.findTracks(sourceLinks, initialParams, ckfOptions);
    if (result.ok()) {
      // Get the track finding output object
      const auto& trackFindingOutput = result.value();
      // Create a SimMultiTrajectory
      trajectories.emplace_back(std::move(trackFindingOutput.fittedStates),
                                std::move(trackFindingOutput.trackTips),
                                std::move(trackFindingOutput.fittedParameters));
    } else {
      ACTS_WARNING("Track finding failed for truth seed "
                   << iseed << " with error" << result.error());
      // Track finding failed, but still create an empty SimMultiTrajectory
      trajectories.push_back(SimMultiTrajectory());
    }
  }

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return ActsExamples::ProcessCode::SUCCESS;
}
