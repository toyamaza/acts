// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <cmath>
#include <vector>

#include <boost/container/static_vector.hpp>

namespace ActsExamples {

/// Space point representation of a measurement suitable for track seeding.
class SimSpacePoint {
 public:
  /// Construct the space point from global position and selected variances.
  ///
  /// @tparam position_t Input position type
  /// @param pos Global position
  /// @param varRho Measurement variance of the global transverse distance
  /// @param varZ Measurement variance of the global longitudinal position
  /// @param measurementIndex Index of the underlying measurement
  template <typename position_t>
  SimSpacePoint(
      const Eigen::MatrixBase<position_t>& pos, float varRho, float varZ,
      boost::container::static_vector<const Acts::SourceLink*, 2> sourceLinks)
      : m_x(pos[Acts::ePos0]),
        m_y(pos[Acts::ePos1]),
        m_z(pos[Acts::ePos2]),
        m_rho(std::hypot(m_x, m_y)),
        m_varianceRho(varRho),
        m_varianceZ(varZ),
        m_sourceLinks(sourceLinks) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(position_t, 3);
  }

  constexpr float x() const { return m_x; }
  constexpr float y() const { return m_y; }
  constexpr float z() const { return m_z; }
  constexpr float r() const { return m_rho; }
  constexpr float varianceR() const { return m_varianceRho; }
  constexpr float varianceZ() const { return m_varianceZ; }

  const boost::container::static_vector<const Acts::SourceLink*, 2>
  sourceLinks() const {
    return m_sourceLinks;
  }

 private:
  // Global position
  float m_x;
  float m_y;
  float m_z;
  float m_rho;
  // Variance in rho/z of the global coordinates
  float m_varianceRho;
  float m_varianceZ;
  // SourceLinks of the corresponding measurements. A Pixel (strip) SP has one
  // (two) sourceLink(s).
  boost::container::static_vector<const Acts::SourceLink*, 2> m_sourceLinks;
};

inline bool operator==(const SimSpacePoint& lhs, const SimSpacePoint& rhs) {
  // TODO would it be sufficient to check just the index under the assumption
  //   that the same measurement index always produces the same space point?
  // no need to check r since it is fully defined by x/y

  return ((lhs.sourceLinks() == rhs.sourceLinks()) and (lhs.x() == rhs.x()) and
          (lhs.y() == rhs.y()) and (lhs.z() == rhs.z()) and
          (lhs.varianceR() == rhs.varianceR()) and
          (lhs.varianceZ() == rhs.varianceZ()));
}

/// Container of space points.
using SimSpacePointContainer = std::vector<SimSpacePoint>;

}  // namespace ActsExamples
