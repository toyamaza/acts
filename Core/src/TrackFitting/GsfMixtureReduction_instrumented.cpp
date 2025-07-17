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
#include <fstream>

// Define this to enable detailed timing instrumentation
#define ENABLE_GSF_TIMING 1

#if ENABLE_GSF_TIMING
class GsfTimer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::string name;
    bool active = false;
    
public:
    GsfTimer(const std::string& timer_name) : name(timer_name) {}
    
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
        active = true;
    }
    
    double stop() {
        if (!active) return 0.0;
        auto end_time = std::chrono::high_resolution_clock::now();
        active = false;
        return std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    }
    
    double elapsed() const {
        if (!active) return 0.0;
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(now - start_time).count();
    }
};

struct TimingResults {
    double matrix_construction = 0.0;
    double minpair_search = 0.0;
    double component_merging = 0.0;
    double distance_recomputation = 0.0;
    double masking = 0.0;
    double cleanup = 0.0;
    double total_time = 0.0;
    int iterations = 0;
    size_t initial_components = 0;
    size_t final_components = 0;
    
    void print() const {
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "=== GSF Mixture Reduction Performance Report ===" << std::endl;
        std::cout << "Components: " << initial_components << " -> " << final_components << " (reduced by " << (initial_components - final_components) << ")" << std::endl;
        std::cout << "Iterations: " << iterations << std::endl;
        std::cout << "Matrix construction:    " << std::setw(8) << matrix_construction << " μs (" << std::setw(5) << (100.0 * matrix_construction / total_time) << "%)" << std::endl;
        std::cout << "Min pair search:        " << std::setw(8) << minpair_search << " μs (" << std::setw(5) << (100.0 * minpair_search / total_time) << "%)" << std::endl;
        std::cout << "Component merging:      " << std::setw(8) << component_merging << " μs (" << std::setw(5) << (100.0 * component_merging / total_time) << "%)" << std::endl;
        std::cout << "Distance recomputation: " << std::setw(8) << distance_recomputation << " μs (" << std::setw(5) << (100.0 * distance_recomputation / total_time) << "%)" << std::endl;
        std::cout << "Masking:                " << std::setw(8) << masking << " μs (" << std::setw(5) << (100.0 * masking / total_time) << "%)" << std::endl;
        std::cout << "Final cleanup:          " << std::setw(8) << cleanup << " μs (" << std::setw(5) << (100.0 * cleanup / total_time) << "%)" << std::endl;
        std::cout << "Total time:             " << std::setw(8) << total_time << " μs" << std::endl;
        if (iterations > 0) {
            std::cout << "Average per iteration:  " << std::setw(8) << ((minpair_search + component_merging + distance_recomputation + masking) / iterations) << " μs" << std::endl;
        }
        std::cout << "=================================================" << std::endl;
    }
    
    void saveToFile(const std::string& filename) const {
        std::ofstream file(filename, std::ios::app);
        if (file.is_open()) {
            file << std::fixed << std::setprecision(3);
            file << initial_components << "," << final_components << "," << iterations << ","
                 << matrix_construction << "," << minpair_search << "," << component_merging << ","
                 << distance_recomputation << "," << masking << "," << cleanup << "," << total_time << std::endl;
        }
    }
};

#define TIMING_START(timer) timer.start()
#define TIMING_STOP(timer) timer.stop()
#define TIMING_ACCUMULATE(var, timer) var += timer.stop()

#else
// No-op macros when timing is disabled
#define TIMING_START(timer)
#define TIMING_STOP(timer) 0.0
#define TIMING_ACCUMULATE(var, timer)
#endif

template <typename proj_t, typename angle_desc_t>
void reduceWithKLDistanceImpl(std::vector<Acts::GsfComponent> &cmpCache,
                              std::size_t maxCmpsAfterMerge, const proj_t &proj,
                              const angle_desc_t &desc) {
#if ENABLE_GSF_TIMING
  TimingResults results;
  results.initial_components = cmpCache.size();
  
  GsfTimer total_timer("total");
  GsfTimer matrix_timer("matrix");
  GsfTimer minpair_timer("minpair");
  GsfTimer merge_timer("merge");
  GsfTimer recompute_timer("recompute");
  GsfTimer mask_timer("mask");
  GsfTimer cleanup_timer("cleanup");
  
  TIMING_START(total_timer);
#endif

  // Timing: Distance matrix construction
#if ENABLE_GSF_TIMING
  TIMING_START(matrix_timer);
#endif
  Acts::detail::SymmetricKLDistanceMatrix distances(cmpCache, proj);
#if ENABLE_GSF_TIMING
  results.matrix_construction = TIMING_STOP(matrix_timer);
#endif

  auto remainingComponents = cmpCache.size();

  while (remainingComponents > maxCmpsAfterMerge) {
#if ENABLE_GSF_TIMING
    results.iterations++;
#endif
    
    // Time finding minimum distance pair
#if ENABLE_GSF_TIMING
    TIMING_START(minpair_timer);
#endif
    const auto [minI, minJ] = distances.minDistancePair();
#if ENABLE_GSF_TIMING
    TIMING_ACCUMULATE(results.minpair_search, minpair_timer);
#endif

    // Time component merging
#if ENABLE_GSF_TIMING
    TIMING_START(merge_timer);
#endif
    cmpCache[minI] =
        mergeComponents(cmpCache[minI], cmpCache[minJ], proj, desc);
#if ENABLE_GSF_TIMING
    TIMING_ACCUMULATE(results.component_merging, merge_timer);
#endif

    // Time distance recomputation
#if ENABLE_GSF_TIMING
    TIMING_START(recompute_timer);
#endif
    distances.recomputeAssociatedDistances(minI, cmpCache, proj);
#if ENABLE_GSF_TIMING
    TIMING_ACCUMULATE(results.distance_recomputation, recompute_timer);
#endif

    // Time masking
#if ENABLE_GSF_TIMING
    TIMING_START(mask_timer);
#endif
    proj(cmpCache[minJ]).weight = -1.0;
    distances.maskAssociatedDistances(minJ);
#if ENABLE_GSF_TIMING
    TIMING_ACCUMULATE(results.masking, mask_timer);
#endif

    remainingComponents--;
  }

  // Timing: Final cleanup
#if ENABLE_GSF_TIMING
  TIMING_START(cleanup_timer);
#endif
  std::ranges::sort(cmpCache, {},
                    [&](const auto &c) { return proj(c).weight; });
  cmpCache.erase(
      std::remove_if(cmpCache.begin(), cmpCache.end(),
                     [&](const auto &a) { return proj(a).weight == -1.0; }),
      cmpCache.end());
#if ENABLE_GSF_TIMING
  results.cleanup = TIMING_STOP(cleanup_timer);
#endif

#if ENABLE_GSF_TIMING
  results.total_time = TIMING_STOP(total_timer);
  results.final_components = cmpCache.size();
  
  // Print and save results
  results.print();
  results.saveToFile("gsf_timing_results.csv");
#endif

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
