// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GsfMixtureReduction.hpp"

#include "Acts/TrackFitting/detail/SymmetricKlDistanceMatrix.hpp"

#include <algorithm>

template <typename proj_t, typename angle_desc_t>
void reduceWithKLDistanceImpl(std::vector<Acts::GsfComponent> &cmpCache,
                              std::size_t maxCmpsAfterMerge, const proj_t &proj,
                              const angle_desc_t &desc) {
  Acts::detail::SymmetricKLDistanceMatrix distances(cmpCache, proj);

  auto remainingComponents = cmpCache.size();

  while (remainingComponents > maxCmpsAfterMerge) {
    auto [minI_current, minJ_current] = distances.minDistancePair();
    
    // Ensure minI is the smaller index to keep indices stable during removal
    if (minI_current > minJ_current) {
        std::swap(minI_current, minJ_current);
    }

    // Get original indices before any potential matrix modification
    const auto minI_original = distances.getOriginalIndex(minI_current);
    const auto minJ_original = distances.getOriginalIndex(minJ_current);

    cmpCache[minI_original] =
        mergeComponents(cmpCache[minI_original], cmpCache[minJ_original], proj, desc);

    // Mark the removed component as invalid for the final cleanup step
    proj(cmpCache[minJ_original]).weight = -1.0;

    distances.recomputeAssociatedDistances(minI_current, cmpCache, proj);

    // Remove the second component from the distance matrix
    // This is the core of the immediate-removal optimization
    distances.removeComponent(minJ_current);

    remainingComponents--;
  }

  // Remove all components which are labeled with weight -1
  std::ranges::sort(cmpCache, {},
                    [&](const auto &c) { return proj(c).weight; });
  cmpCache.erase(
      std::remove_if(cmpCache.begin(), cmpCache.end(),
                     [&](const auto &a) { return proj(a).weight == -1.0; }),
      cmpCache.end());

  assert(cmpCache.size() == maxCmpsAfterMerge && "size mismatch");
}

namespace Acts {

void reduceMixtureLargestWeights(std::vector<GsfComponent> &cmpCache,
                                 std::size_t maxCmpsAfterMerge,
                                 const Surface & /*surface*/) {
  if (cmpCache.size() <= maxCmpsAfterMerge) {
    return;
  }

  std::nth_element(
      cmpCache.begin(), cmpCache.begin() + maxCmpsAfterMerge, cmpCache.end(),
      [](const auto &a, const auto &b) { return a.weight > b.weight; });
  cmpCache.resize(maxCmpsAfterMerge);
}

void reduceMixtureWithKLDistance(std::vector<Acts::GsfComponent> &cmpCache,
                                 std::size_t maxCmpsAfterMerge,
                                 const Surface &surface) {
  if (cmpCache.size() <= maxCmpsAfterMerge) {
    return;
  }

  auto proj = [](auto &a) -> decltype(auto) { return a; };

  // We must differ between surface types, since there can be different
  // local coordinates
  detail::angleDescriptionSwitch(surface, [&](const auto &desc) {
    reduceWithKLDistanceImpl(cmpCache, maxCmpsAfterMerge, proj, desc);
  });
}

}  // namespace Acts
