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
#include <chrono>
#include <iostream>
#include <iomanip>
#include <mutex>
#include <vector>

// Global timing statistics collector
struct GsfTimingStats {
    std::mutex mutex;
    std::vector<double> matrix_construction_times;
    std::vector<double> minpair_times;
    std::vector<double> merge_times;
    std::vector<double> recompute_times;
    std::vector<double> mask_times;
    std::vector<double> cleanup_times;
    std::vector<double> total_times;
    std::vector<size_t> initial_components;
    std::vector<size_t> final_components;
    std::vector<int> iterations;
    
    void addTiming(double matrix_time, double minpair_time, double merge_time,
                   double recompute_time, double mask_time, double cleanup_time,
                   double total_time, size_t initial_comp, size_t final_comp, int iter) {
        std::lock_guard<std::mutex> lock(mutex);
        matrix_construction_times.push_back(matrix_time);
        minpair_times.push_back(minpair_time);
        merge_times.push_back(merge_time);
        recompute_times.push_back(recompute_time);
        mask_times.push_back(mask_time);
        cleanup_times.push_back(cleanup_time);
        total_times.push_back(total_time);
        initial_components.push_back(initial_comp);
        final_components.push_back(final_comp);
        iterations.push_back(iter);
    }
    
    void printSummary() {
        std::lock_guard<std::mutex> lock(mutex);
        if (total_times.empty()) {
            std::cout << "No timing data collected yet." << std::endl;
            return;
        }
        
        size_t num_calls = total_times.size();
        
        auto sum = [](const std::vector<double>& v) {
            double total = 0.0;
            for (double val : v) total += val;
            return total;
        };
        
        auto avg = [&](const std::vector<double>& v) {
            return sum(v) / v.size();
        };
        
        auto sum_int = [](const std::vector<int>& v) {
            int total = 0;
            for (int val : v) total += val;
            return total;
        };
        
        auto sum_size = [](const std::vector<size_t>& v) {
            size_t total = 0;
            for (size_t val : v) total += val;
            return total;
        };
        
        double total_sum = sum(total_times);
        double matrix_sum = sum(matrix_construction_times);
        double minpair_sum = sum(minpair_times);
        double merge_sum = sum(merge_times);
        double recompute_sum = sum(recompute_times);
        double mask_sum = sum(mask_times);
        double cleanup_sum = sum(cleanup_times);
        
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "\n=== GSF Mixture Reduction AGGREGATE Performance Report ===" << std::endl;
        std::cout << "Total function calls: " << num_calls << std::endl;
        std::cout << "Total components processed: " << sum_size(initial_components) << " -> " << sum_size(final_components) << std::endl;
        std::cout << "Total iterations: " << sum_int(iterations) << std::endl;
        std::cout << "Average components per call: " << (sum_size(initial_components) / num_calls) << " -> " << (sum_size(final_components) / num_calls) << std::endl;
        std::cout << "Average iterations per call: " << (sum_int(iterations) / double(num_calls)) << std::endl;
        std::cout << std::endl;
        
        std::cout << "CUMULATIVE TIMES (across all calls):" << std::endl;
        std::cout << "Matrix construction:    " << std::setw(10) << matrix_sum << " μs (" << std::setw(5) << (100.0 * matrix_sum / total_sum) << "%)" << std::endl;
        std::cout << "Min pair search:        " << std::setw(10) << minpair_sum << " μs (" << std::setw(5) << (100.0 * minpair_sum / total_sum) << "%)" << std::endl;
        std::cout << "Component merging:      " << std::setw(10) << merge_sum << " μs (" << std::setw(5) << (100.0 * merge_sum / total_sum) << "%)" << std::endl;
        std::cout << "Distance recomputation: " << std::setw(10) << recompute_sum << " μs (" << std::setw(5) << (100.0 * recompute_sum / total_sum) << "%)" << std::endl;
        std::cout << "Masking:                " << std::setw(10) << mask_sum << " μs (" << std::setw(5) << (100.0 * mask_sum / total_sum) << "%)" << std::endl;
        std::cout << "Final cleanup:          " << std::setw(10) << cleanup_sum << " μs (" << std::setw(5) << (100.0 * cleanup_sum / total_sum) << "%)" << std::endl;
        std::cout << "TOTAL TIME:             " << std::setw(10) << total_sum << " μs" << std::endl;
        std::cout << std::endl;
        
        std::cout << "AVERAGE TIMES (per call):" << std::endl;
        std::cout << "Matrix construction:    " << std::setw(10) << avg(matrix_construction_times) << " μs" << std::endl;
        std::cout << "Min pair search:        " << std::setw(10) << avg(minpair_times) << " μs" << std::endl;
        std::cout << "Component merging:      " << std::setw(10) << avg(merge_times) << " μs" << std::endl;
        std::cout << "Distance recomputation: " << std::setw(10) << avg(recompute_times) << " μs" << std::endl;
        std::cout << "Masking:                " << std::setw(10) << avg(mask_times) << " μs" << std::endl;
        std::cout << "Final cleanup:          " << std::setw(10) << avg(cleanup_times) << " μs" << std::endl;
        std::cout << "TOTAL TIME:             " << std::setw(10) << avg(total_times) << " μs" << std::endl;
        std::cout << std::endl;
        
        // Performance insights
        std::cout << "PERFORMANCE INSIGHTS:" << std::endl;
        if (recompute_sum > 0.4 * total_sum) {
            std::cout << "⚠️  Distance recomputation is the major bottleneck (" << (100.0 * recompute_sum / total_sum) << "%)" << std::endl;
        }
        if (minpair_sum > 0.3 * total_sum) {
            std::cout << "⚠️  Min pair search is expensive (" << (100.0 * minpair_sum / total_sum) << "%)" << std::endl;
        }
        if (matrix_sum > 0.2 * total_sum) {
            std::cout << "⚠️  Matrix construction overhead is significant (" << (100.0 * matrix_sum / total_sum) << "%)" << std::endl;
        }
        std::cout << "=============================================================" << std::endl;
    }
};

// Global instance
static GsfTimingStats g_timing_stats;

template <typename proj_t, typename angle_desc_t>
void reduceWithKLDistanceImpl(std::vector<Acts::GsfComponent> &cmpCache,
                              std::size_t maxCmpsAfterMerge, const proj_t &proj,
                              const angle_desc_t &desc) {
  auto start_total = std::chrono::high_resolution_clock::now();
  size_t initial_size = cmpCache.size();
  
  // Timing: Distance matrix construction
  auto start_matrix = std::chrono::high_resolution_clock::now();
  Acts::detail::SymmetricKLDistanceMatrix distances(cmpCache, proj);
  auto end_matrix = std::chrono::high_resolution_clock::now();
  auto matrix_time = std::chrono::duration_cast<std::chrono::microseconds>(end_matrix - start_matrix).count();

  auto remainingComponents = cmpCache.size();

  // Timing: Main reduction loop
  double minpair_time = 0.0;
  double merge_time = 0.0;
  double recompute_time = 0.0;
  double mask_time = 0.0;
  int iterations = 0;

  while (remainingComponents > maxCmpsAfterMerge) {
    iterations++;
    
    // Time finding minimum distance pair
    auto start_minpair = std::chrono::high_resolution_clock::now();
    const auto [minI, minJ] = distances.minDistancePair();
    auto end_minpair = std::chrono::high_resolution_clock::now();
    minpair_time += std::chrono::duration_cast<std::chrono::microseconds>(end_minpair - start_minpair).count();

    // Time component merging
    auto start_merge = std::chrono::high_resolution_clock::now();
    cmpCache[minI] =
        mergeComponents(cmpCache[minI], cmpCache[minJ], proj, desc);
    auto end_merge = std::chrono::high_resolution_clock::now();
    merge_time += std::chrono::duration_cast<std::chrono::microseconds>(end_merge - start_merge).count();

    // Time distance recomputation
    auto start_recompute = std::chrono::high_resolution_clock::now();
    distances.recomputeAssociatedDistances(minI, cmpCache, proj);
    auto end_recompute = std::chrono::high_resolution_clock::now();
    recompute_time += std::chrono::duration_cast<std::chrono::microseconds>(end_recompute - start_recompute).count();

    // Time masking
    auto start_mask = std::chrono::high_resolution_clock::now();
    proj(cmpCache[minJ]).weight = -1.0;
    distances.maskAssociatedDistances(minJ);
    auto end_mask = std::chrono::high_resolution_clock::now();
    mask_time += std::chrono::duration_cast<std::chrono::microseconds>(end_mask - start_mask).count();

    remainingComponents--;
  }

  // Timing: Final cleanup - optimized with partition instead of sort
  auto start_cleanup = std::chrono::high_resolution_clock::now();
  // Use partition for O(n) instead of sort for O(n log n)
  auto partition_point = std::partition(cmpCache.begin(), cmpCache.end(),
                                       [&](const auto &a) { return proj(a).weight != -1.0; });
  cmpCache.erase(partition_point, cmpCache.end());
  auto end_cleanup = std::chrono::high_resolution_clock::now();
  auto cleanup_time = std::chrono::duration_cast<std::chrono::microseconds>(end_cleanup - start_cleanup).count();

  auto end_total = std::chrono::high_resolution_clock::now();
  auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end_total - start_total).count();

  // Store timing data for aggregate analysis
  g_timing_stats.addTiming(matrix_time, minpair_time, merge_time, recompute_time, 
                          mask_time, cleanup_time, total_time, 
                          initial_size, cmpCache.size(), iterations);

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

// Function to print aggregate timing statistics
void printGsfTimingSummary() {
  g_timing_stats.printSummary();
}

// Function to reset timing statistics
void resetGsfTimingStats() {
  std::lock_guard<std::mutex> lock(g_timing_stats.mutex);
  g_timing_stats.matrix_construction_times.clear();
  g_timing_stats.minpair_times.clear();
  g_timing_stats.merge_times.clear();
  g_timing_stats.recompute_times.clear();
  g_timing_stats.mask_times.clear();
  g_timing_stats.cleanup_times.clear();
  g_timing_stats.total_times.clear();
  g_timing_stats.initial_components.clear();
  g_timing_stats.final_components.clear();
  g_timing_stats.iterations.clear();
}

}  // namespace Acts
