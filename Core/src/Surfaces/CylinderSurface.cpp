// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/CylinderSurface.hpp"

#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <cassert>
#include <cmath>
#include <system_error>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::CylinderSurface::CylinderSurface(const CylinderSurface& other)
    : GeometryObject(), Surface(other), m_bounds(other.m_bounds) {}

Acts::CylinderSurface::CylinderSurface(const GeometryContext& gctx,
                                       const CylinderSurface& other,
                                       const Transform3D& shift)
    : GeometryObject(), Surface(gctx, other, shift), m_bounds(other.m_bounds) {}

Acts::CylinderSurface::CylinderSurface(const Transform3D& transform,
                                       double radius, double halfz,
                                       double halfphi, double avphi)
    : Surface(transform),
      m_bounds(std::make_shared<const CylinderBounds>(radius, halfz, halfphi,
                                                      avphi)) {}

Acts::CylinderSurface::CylinderSurface(
    std::shared_ptr<const CylinderBounds> cbounds,
    const Acts::DetectorElementBase& detelement)
    : Surface(detelement), m_bounds(std::move(cbounds)) {
  /// surfaces representing a detector element must have bounds
  assert(cbounds);
}

Acts::CylinderSurface::CylinderSurface(
    const Transform3D& transform,
    const std::shared_ptr<const CylinderBounds>& cbounds)
    : Surface(transform), m_bounds(cbounds) {
  throw_assert(cbounds, "CylinderBounds must not be nullptr");
}

Acts::CylinderSurface& Acts::CylinderSurface::operator=(
    const CylinderSurface& other) {
  if (this != &other) {
    Surface::operator=(other);
    m_bounds = other.m_bounds;
  }
  return *this;
}

// return the binning position for ordering in the BinnedArray
Acts::Vector3D Acts::CylinderSurface::binningPosition(
    const GeometryContext& gctx, BinningValue bValue) const {
  const Acts::Vector3D& sfCenter = center(gctx);
  // special binning type for R-type methods
  if (bValue == Acts::binR || bValue == Acts::binRPhi) {
    double R = bounds().get(CylinderBounds::eR);
    double phi = bounds().get(CylinderBounds::eAveragePhi);
    return Vector3D(sfCenter.x() + R * cos(phi), sfCenter.y() + R * sin(phi),
                    sfCenter.z());
  }
  // give the center as default for all of these binning types
  // binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta
  return sfCenter;
}

// return the measurement frame: it's the tangential plane
Acts::RotationMatrix3D Acts::CylinderSurface::referenceFrame(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& /*unused*/) const {
  RotationMatrix3D mFrame;
  // construct the measurement frame
  // measured Y is the z axis
  Vector3D measY = rotSymmetryAxis(gctx);
  // measured z is the position normalized transverse (in local)
  Vector3D measDepth = normal(gctx, position);
  // measured X is what comoes out of it
  Vector3D measX(measY.cross(measDepth).normalized());
  // assign the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  return mFrame;
}

Acts::Surface::SurfaceType Acts::CylinderSurface::type() const {
  return Surface::Cylinder;
}

Acts::Vector3D Acts::CylinderSurface::localToGlobal(
    const GeometryContext& gctx, const Vector2D& lposition,
    const Vector3D& /*unused*/) const {
  // create the position in the local 3d frame
  double r = bounds().get(CylinderBounds::eR);
  double phi = lposition[Acts::eBoundLoc0] / r;
  Vector3D position(r * cos(phi), r * sin(phi), lposition[Acts::eBoundLoc1]);
  return transform(gctx) * position;
}

Acts::Result<Acts::Vector2D> Acts::CylinderSurface::globalToLocal(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& /*unused*/) const {
  // @todo check if s_onSurfaceTolerance would do here
  double inttol = bounds().get(CylinderBounds::eR) * 0.0001;
  if (inttol < 0.01) {
    inttol = 0.01;
  }
  const Transform3D& sfTransform = transform(gctx);
  Transform3D inverseTrans(sfTransform.inverse());
  Vector3D loc3Dframe(inverseTrans * position);
  if (std::abs(perp(loc3Dframe) - bounds().get(CylinderBounds::eR)) > inttol) {
    return Result<Vector2D>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }
  return Result<Vector2D>::success(
      {bounds().get(CylinderBounds::eR) * phi(loc3Dframe), loc3Dframe.z()});
}

std::string Acts::CylinderSurface::name() const {
  return "Acts::CylinderSurface";
}

Acts::Vector3D Acts::CylinderSurface::normal(
    const GeometryContext& gctx, const Acts::Vector2D& lposition) const {
  double phi = lposition[Acts::eBoundLoc0] / m_bounds->get(CylinderBounds::eR);
  Vector3D localNormal(cos(phi), sin(phi), 0.);
  return Vector3D(transform(gctx).matrix().block<3, 3>(0, 0) * localNormal);
}

Acts::Vector3D Acts::CylinderSurface::normal(
    const GeometryContext& gctx, const Acts::Vector3D& position) const {
  const Transform3D& sfTransform = transform(gctx);
  // get it into the cylinder frame
  Vector3D pos3D = sfTransform.inverse() * position;
  // set the z coordinate to 0
  pos3D.z() = 0.;
  // normalize and rotate back into global if needed
  return sfTransform.linear() * pos3D.normalized();
}

double Acts::CylinderSurface::pathCorrection(
    const GeometryContext& gctx, const Acts::Vector3D& position,
    const Acts::Vector3D& direction) const {
  Vector3D normalT = normal(gctx, position);
  double cosAlpha = normalT.dot(direction);
  return std::fabs(1. / cosAlpha);
}

const Acts::CylinderBounds& Acts::CylinderSurface::bounds() const {
  return (*m_bounds.get());
}

Acts::Polyhedron Acts::CylinderSurface::polyhedronRepresentation(
    const GeometryContext& gctx, size_t lseg) const {
  // Prepare vertices and faces
  std::vector<Vector3D> vertices;
  std::vector<Polyhedron::FaceType> faces;
  std::vector<Polyhedron::FaceType> triangularMesh;

  auto ctrans = transform(gctx);
  bool fullCylinder = bounds().coversFullAzimuth();

  double avgPhi = bounds().get(CylinderBounds::eAveragePhi);
  double halfPhi = bounds().get(CylinderBounds::eHalfPhiSector);

  // Get the phi segments from the helper - ensures extra points
  auto phiSegs = fullCylinder
                     ? detail::VerticesHelper::phiSegments()
                     : detail::VerticesHelper::phiSegments(
                           avgPhi - halfPhi, avgPhi + halfPhi, {avgPhi});

  // Write the two bows/circles on either side
  std::vector<int> sides = {-1, 1};
  for (auto& side : sides) {
    for (size_t iseg = 0; iseg < phiSegs.size() - 1; ++iseg) {
      int addon = (iseg == phiSegs.size() - 2 and not fullCylinder) ? 1 : 0;
      /// Helper method to create the segment
      detail::VerticesHelper::createSegment(
          vertices,
          {bounds().get(CylinderBounds::eR), bounds().get(CylinderBounds::eR)},
          phiSegs[iseg], phiSegs[iseg + 1], lseg, addon,
          Vector3D(0., 0., side * bounds().get(CylinderBounds::eHalfLengthZ)),
          ctrans);
    }
  }
  auto facesMesh =
      detail::FacesHelper::cylindricalFaceMesh(vertices, fullCylinder);
  return Polyhedron(vertices, facesMesh.first, facesMesh.second, false);
}
