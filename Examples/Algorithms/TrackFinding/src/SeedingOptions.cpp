// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingOptions.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <string>

#include <boost/program_options.hpp>

void ActsExamples::Options::addSeedingOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("seed-config-file", value<std::string>()->default_value(""),
      "Configuration (.json) file for seeding");
  
}

ActsExamples::SeedingAlgorithm::Config ActsExamples::Options::readSeedingConfig(
    const ActsExamples::Options::Variables& variables) {
  auto chi2Max = variables["ckf-selection-chi2max"].template as<double>();
  auto nMax = variables["ckf-selection-nmax"].template as<size_t>();

  // config is a GeometryHierarchyMap with just the global default
  SeedingAlgorithm::Config cfg;

  return cfg;
}
