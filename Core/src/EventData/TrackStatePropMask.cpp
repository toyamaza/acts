// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackStatePropMask.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <ostream>

namespace Acts {

std::ostream& operator<<(std::ostream& os, TrackStatePropMask mask) {
  using PM = TrackStatePropMask;
  os << "TrackStatePropMask(";
  if (mask == PM::None) {
    os << "None";
  } else {
    os << "\n  [" << (ACTS_CHECK_BIT(mask, PM::Predicted) ? "x" : " ")
       << "] predicted";
    os << "\n  [" << (ACTS_CHECK_BIT(mask, PM::Filtered) ? "x" : " ")
       << "] filtered";
    os << "\n  [" << (ACTS_CHECK_BIT(mask, PM::Smoothed) ? "x" : " ")
       << "] smoothed";
    os << "\n  [" << (ACTS_CHECK_BIT(mask, PM::Jacobian) ? "x" : " ")
       << "] jacobian";
    os << "\n  [" << (ACTS_CHECK_BIT(mask, PM::Uncalibrated) ? "x" : " ")
       << "] uncalibrated";
    os << "\n  [" << (ACTS_CHECK_BIT(mask, PM::Calibrated) ? "x" : " ")
       << "] calibrated";
    os << "\n";
  }
  os << ")";
  return os;
}

}  // namespace Acts
