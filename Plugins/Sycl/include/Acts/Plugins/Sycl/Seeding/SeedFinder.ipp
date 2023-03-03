// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// System include(s)
#include <algorithm>
#include <cmath>
#include <utility>

// VecMem include(s).
#include "vecmem/containers/vector.hpp"

// Acts include(s).
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"

// SYCL plugin include(s)
#include "Acts/Plugins/Sycl/Seeding/CreateSeedsForGroupSycl.hpp"
#include "Acts/Plugins/Sycl/Seeding/SeedFinder.hpp"

namespace Acts::Sycl {
template <typename external_spacepoint_t>
SeedFinder<external_spacepoint_t>::SeedFinder(
    Acts::SeedFinderConfig<external_spacepoint_t> config,
    const Acts::SeedFinderOptions& options,
    const Acts::Sycl::DeviceExperimentCuts& cuts,
    Acts::Sycl::QueueWrapper wrappedQueue, vecmem::memory_resource& resource,
    vecmem::memory_resource* device_resource)
    : m_config(config),
      m_options(options),
      m_deviceCuts(cuts),
      m_wrappedQueue(std::move(wrappedQueue)),
      m_resource(&resource),
      m_device_resource(device_resource) {
  auto seedFilterConfig = m_config.seedFilter->getSeedFilterConfig();

  // init m_deviceConfig
  m_deviceConfig = Acts::Sycl::detail::DeviceSeedFinderConfig{
      m_config.deltaRMin,
      m_config.deltaRMax,
      m_config.cotThetaMax,
      m_config.collisionRegionMin,
      m_config.collisionRegionMax,
      m_config.maxScatteringAngle2,
      m_config.sigmaScattering,
      m_options.minHelixDiameter2,
      m_options.pT2perRadius,
      seedFilterConfig.deltaInvHelixDiameter,
      seedFilterConfig.impactWeightFactor,
      seedFilterConfig.deltaRMin,
      seedFilterConfig.compatSeedWeight,
      m_config.impactMax,
      seedFilterConfig.compatSeedLimit,
  };
}

template <typename external_spacepoint_t>
template <typename sp_range_t>
std::vector<Acts::Seed<external_spacepoint_t>>
SeedFinder<external_spacepoint_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  std::vector<Seed<external_spacepoint_t>> outputVec;

  // As a first step, we create Arrays of Structures (AoS)
  // that are easily comprehensible by the GPU. This allows us
  // less memory access operations than with simple (float) arrays.

  // Creating VecMem vectors of the space points, linked to the host/shared
  // memory resource They will be filled and passed to CreateSeedsForGroup().
  vecmem::vector<detail::DeviceSpacePoint> deviceBottomSPs(m_resource);
  vecmem::vector<detail::DeviceSpacePoint> deviceMiddleSPs(m_resource);
  vecmem::vector<detail::DeviceSpacePoint> deviceTopSPs(m_resource);

  std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*> bottomSPvec;
  std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*> middleSPvec;
  std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*> topSPvec;

  for (auto SP : bottomSPs) {
    bottomSPvec.push_back(SP);
  }
  deviceBottomSPs.reserve(bottomSPvec.size());
  for (auto SP : bottomSPvec) {
    deviceBottomSPs.push_back({SP->x(), SP->y(), SP->z(), SP->radius(),
                               SP->varianceR(), SP->varianceZ()});
  }

  for (auto SP : middleSPs) {
    middleSPvec.push_back(SP);
  }
  deviceMiddleSPs.reserve(middleSPvec.size());
  for (auto SP : middleSPvec) {
    deviceMiddleSPs.push_back({SP->x(), SP->y(), SP->z(), SP->radius(),
                               SP->varianceR(), SP->varianceZ()});
  }

  for (auto SP : topSPs) {
    topSPvec.push_back(SP);
  }
  deviceTopSPs.reserve(topSPvec.size());
  for (auto SP : topSPvec) {
    deviceTopSPs.push_back({SP->x(), SP->y(), SP->z(), SP->radius(),
                            SP->varianceR(), SP->varianceZ()});
  }

  // std::vector<std::vector<detail::SeedData>> seeds;
  std::vector<std::vector<detail::SeedData>> seeds;

  // Call the SYCL seeding algorithm
  createSeedsForGroupSycl(m_wrappedQueue, *m_resource, m_device_resource,
                          m_deviceConfig, m_deviceCuts, deviceBottomSPs,
                          deviceMiddleSPs, deviceTopSPs, seeds);

  // Iterate through seeds returned by the SYCL algorithm and perform the last
  // step of filtering for fixed middle SP.
  std::vector<std::pair<
      float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
      seedsPerSPM;
  for (size_t mi = 0; mi < seeds.size(); ++mi) {
    seedsPerSPM.clear();
    for (size_t j = 0; j < seeds[mi].size(); ++j) {
      auto& bottomSP = *(bottomSPvec[seeds[mi][j].bottom]);
      auto& middleSP = *(middleSPvec[mi]);
      auto& topSP = *(topSPvec[seeds[mi][j].top]);
      float weight = seeds[mi][j].weight;

      seedsPerSPM.emplace_back(std::make_pair(
          weight, std::make_unique<const InternalSeed<external_spacepoint_t>>(
                      bottomSP, middleSP, topSP, 0)));
    }
    int numQualitySeeds = 0;  // not used but needs to be fixed
    m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSPM, numQualitySeeds,
                                              std::back_inserter(outputVec));
  }
  return outputVec;
}
}  // namespace Acts::Sycl
