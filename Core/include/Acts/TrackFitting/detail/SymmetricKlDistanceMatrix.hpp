// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <numeric>

#include "Acts/EventData/MultiComponentTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/TrackFitting/detail/GsfComponentMerging.hpp"
#include "Acts/TrackFitting/detail/GsfUtils.hpp"

namespace Acts::detail {

/// Computes the Kullback-Leibler distance between two components as shown in
/// https://arxiv.org/abs/2001.00727v1 but ignoring the weights
template <typename component_t, typename component_projector_t>
auto computeSymmetricKlDivergence(const component_t &a, const component_t &b,
                                  const component_projector_t &proj) {
  using namespace Acts;
  const auto parsA = proj(a).boundPars[eBoundQOverP];
  const auto parsB = proj(b).boundPars[eBoundQOverP];
  const auto covA = proj(a).boundCov(eBoundQOverP, eBoundQOverP);
  const auto covB = proj(b).boundCov(eBoundQOverP, eBoundQOverP);

  const auto kl = covA * (1 / covB) + covB * (1 / covA) +
                  (parsA - parsB) * (1 / covA + 1 / covB) * (parsA - parsB);

  return kl;
}

template <typename component_t, typename component_projector_t,
          typename angle_desc_t>
auto mergeComponents(const component_t &a, const component_t &b,
                     const component_projector_t &proj,
                     const angle_desc_t &angle_desc) {
  assert(proj(a).weight >= 0.0 && proj(b).weight >= 0.0 &&
         "non-positive weight");

  std::array range = {std::ref(proj(a)), std::ref(proj(b))};
  const auto refProj = [](auto &c) {
    return std::tie(c.get().weight, c.get().boundPars, c.get().boundCov);
  };

  auto [mergedPars, mergedCov] =
      gaussianMixtureMeanCov(range, refProj, angle_desc);

  component_t ret = a;
  proj(ret).boundPars = mergedPars;
  proj(ret).boundCov = mergedCov;
  proj(ret).weight = proj(a).weight + proj(b).weight;

  return ret;
}

/// @brief Class representing a symmetric distance matrix
class SymmetricKLDistanceMatrix {
  using Array = Eigen::Array<double, Eigen::Dynamic, 1>;

  Array m_distances;
  std::vector<std::pair<std::size_t, std::size_t>> m_mapToPair;
  std::size_t m_numberComponents;
  
  // Mapping from current matrix indices to original component indices
  std::vector<std::size_t> m_currentToOriginal;

 public:
  template <typename component_t, typename projector_t>
  SymmetricKLDistanceMatrix(const std::vector<component_t> &cmps,
                            const projector_t &proj)
      : m_distances(Array::Zero(cmps.size() * (cmps.size() - 1) / 2)),
        m_mapToPair(m_distances.size()),
        m_numberComponents(cmps.size()) {
    // Initialize mapping from current to original indices
    m_currentToOriginal.resize(m_numberComponents);
    std::iota(m_currentToOriginal.begin(), m_currentToOriginal.end(), 0);
    
    for (auto i = 1ul; i < m_numberComponents; ++i) {
      const auto indexConst = (i - 1) * i / 2;
      for (auto j = 0ul; j < i; ++j) {
        m_mapToPair.at(indexConst + j) = {i, j};
        m_distances[indexConst + j] =
            computeSymmetricKlDivergence(cmps[i], cmps[j], proj);
      }
    }
  }

  auto at(std::size_t i, std::size_t j) const {
    return m_distances[i * (i - 1) / 2 + j];
  }

  template <typename component_t, typename projector_t>
  void recomputeAssociatedDistances(std::size_t n,
                                    const std::vector<component_t> &cmps,
                                    const projector_t &proj) {
    const auto original_n = m_currentToOriginal[n];
    
    // Column n (i > n)
    for (std::size_t i = n + 1; i < m_numberComponents; i++) {
        const auto original_i = m_currentToOriginal[i];
        m_distances[(i * (i - 1) / 2) + n] = computeSymmetricKlDivergence(cmps[original_i], cmps[original_n], proj);
    }
    // Row n (j < n)
    for (std::size_t j = 0; j < n; j++) {
        const auto original_j = m_currentToOriginal[j];
        m_distances[(n * (n - 1) / 2) + j] = computeSymmetricKlDivergence(cmps[original_n], cmps[original_j], proj);
    }
  }

  // Remove a component from the matrix by swapping with the last and shrinking (O(N) operation)
  void removeComponent(std::size_t componentToRemove) {
    if (m_numberComponents <= 1 || componentToRemove >= m_numberComponents) return;

    std::size_t lastComponent = m_numberComponents - 1;

    if (componentToRemove != lastComponent) {
        // Copy distances from the last component over the one to be removed
        
        // Copy row `lastComponent` to row `componentToRemove`
        // This handles pairs (componentToRemove, j) where j < componentToRemove
        for (std::size_t j = 0; j < componentToRemove; ++j) {
            m_distances[(componentToRemove * (componentToRemove - 1) / 2) + j] = m_distances[(lastComponent * (lastComponent - 1) / 2) + j];
        }
        
        // Copy column `lastComponent` to column `componentToRemove`
        // This handles pairs (i, componentToRemove) where i > componentToRemove
        for (std::size_t i = componentToRemove + 1; i < lastComponent; ++i) {
            m_distances[(i * (i - 1) / 2) + componentToRemove] = m_distances[(i * (i - 1) / 2) + lastComponent];
        }
        
        // Update the original index mapping
        m_currentToOriginal[componentToRemove] = m_currentToOriginal[lastComponent];
    }
    
    // Logically shrink the matrix
    m_numberComponents--;
    m_currentToOriginal.pop_back();
  }

  // Get the original component index for a current matrix index
  std::size_t getOriginalIndex(std::size_t currentIndex) const {
    return m_currentToOriginal[currentIndex];
  }

  // Return the number of active components
  std::size_t getActiveComponentCount() const {
    return m_numberComponents;
  }

  auto minDistancePair() const {
    auto min = std::numeric_limits<double>::max();
    std::size_t idx = 0;
    
    std::size_t current_dist_size = m_numberComponents * (m_numberComponents - 1) / 2;
    for (auto i = 0l; i < current_dist_size; ++i) {
      if (auto distance = m_distances[i]; distance < min) {
        min = distance;
        idx = i;
      }
    }

    return m_mapToPair.at(idx);
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const SymmetricKLDistanceMatrix &m) {
    return os;
  }
};

}  // namespace Acts::detail
