// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>        // for string printing
#include <system_error>  // bring in std::error_code et al

namespace Acts {
// This is the custom error code enum
enum class SurfaceError { GlobalPositionNotOnSurface = 1 };

namespace detail {
// Define a custom error code category derived from std::error_category
class SurfaceErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category
  const char* name() const noexcept final { return "SurfaceError"; }
  // Return what each enum means in text
  std::string message(int c) const final {
    switch (static_cast<SurfaceError>(c)) {
      case SurfaceError::GlobalPositionNotOnSurface:
        return "Global to local transformation failed: position not on "
               "surface.";
      default:
        return "unknown";
    }
  }
};
}  // namespace detail

// Declare a global function returning a static instance of the custom category
extern inline const detail::SurfaceErrorCategory& SurfaceErrorCategory() {
  static detail::SurfaceErrorCategory c;
  return c;
}

inline std::error_code make_error_code(Acts::SurfaceError e) {
  return {static_cast<int>(e), Acts::SurfaceErrorCategory()};
}
}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::SurfaceError> : std::true_type {};
}  // namespace std