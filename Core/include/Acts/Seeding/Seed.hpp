// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

namespace Acts {
template <typename SpacePoint>
class Seed {
  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

 public:
  Seed(const SpacePoint& b, const SpacePoint& m, const SpacePoint& u,
       float vertex);
  Seed(const Seed&) = default;
  Seed& operator=(const Seed&);

  const std::vector<const SpacePoint*>& sp() const { return m_spacepoints; }
  double z() const { return m_zvertex; }
  double invHelixDiameter() const { return m_invHelixDiameter; }
  double impactParameter() const { return m_impactParameter; }
  double cotTheta() const { return m_cotTheta; }

 private:
  std::vector<const SpacePoint*> m_spacepoints;
  float m_zvertex;
  float m_invHelixDiameter;
  float m_impactParameter;
  float m_cotTheta;
};

///////////////////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
Seed<SpacePoint>::Seed(const SpacePoint& b, const SpacePoint& m,
                       const SpacePoint& u, float vertex) {
  m_zvertex = vertex;
  m_spacepoints.push_back(&b);
  m_spacepoints.push_back(&m);
  m_spacepoints.push_back(&u);
  m_invHelixDiameter = -999.; // Fix me
  m_impactParameter =  -999.; // Fix me
  m_cotTheta =  -999.; // Fix me
}

template <typename SpacePoint>
Seed<SpacePoint>::Seed(const SpacePoint& b, const SpacePoint& m,
                       const SpacePoint& u, float vertex,
		       float invHelixDiameter, float impactParameter,
		       float cotTheta) {
  m_zvertex = vertex;
  m_spacepoints.push_back(&b);
  m_spacepoints.push_back(&m);
  m_spacepoints.push_back(&u);
  m_invHelixDiameter = invHelixDiameter;
  m_impactParameter = impactParameter;
  m_cotTheta = cotTheta;
}



}  // namespace Acts
