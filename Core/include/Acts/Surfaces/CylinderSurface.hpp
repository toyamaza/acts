// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryStatics.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/RealQuadraticEquation.hpp"

#include <cmath>

namespace Acts {

/// @class CylinderSurface
///
/// Class for a CylinderSurface in the TrackingGeometry.
/// It inherits from Surface.
///
/// The cylinder surface has a special role in the TrackingGeometry,
/// since it builds the surfaces of all TrackingVolumes at container level
/// for a cylindrical tracking geometry.
///
/// @image html CylinderSurface.png

class CylinderSurface : public Surface {
  friend Surface;

 protected:
  /// Constructor from DetectorElementBase: Element proxy
  ///
  /// @param cbounds are the provided cylinder bounds (shared)
  /// @param detelement is the linked detector element to this surface
  CylinderSurface(std::shared_ptr<const CylinderBounds> cbounds,
                  const DetectorElementBase& detelement);

  /// Constructor from Transform3D and CylinderBounds
  ///
  /// @param transform The transform to position the surface
  /// @param radius The radius of the cylinder
  /// @param halfz The half length in z
  /// @param halfphi The half opening angle
  /// @param avphi The phi value from which the opening angle spans (both sides)
  CylinderSurface(const Transform3D& transform, double radius, double halfz,
                  double halfphi = M_PI, double avphi = 0.);

  /// Constructor from Transform3D and CylinderBounds arguments
  ///
  /// @param transform The transform to position the surface
  /// @param cbounds is a shared pointer to a cylindeer bounds object,
  /// it must exist (assert test)
  CylinderSurface(const Transform3D& transform,
                  const std::shared_ptr<const CylinderBounds>& cbounds);

  /// Copy constructor
  ///
  /// @param other is the source cylinder for the copy
  CylinderSurface(const CylinderSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other is the source cone surface
  /// @param shift is the additional transfrom applied after copying
  CylinderSurface(const GeometryContext& gctx, const CylinderSurface& other,
                  const Transform3D& shift);

 public:
  ~CylinderSurface() override = default;
  CylinderSurface() = delete;

  /// Assignment operator
  ///
  /// @param other is the source cylinder for the copy
  CylinderSurface& operator=(const CylinderSurface& other);

  /// The binning position method - is overloaded for r-type binning
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the type of global binning to be done
  ///
  /// @return is the global position to be used for binning
  Vector3D binningPosition(const GeometryContext& gctx,
                           BinningValue bValue) const final;

  /// Return the measurement frame - this is needed for alignment, in particular
  /// The measurement frame of a cylinder is the tangential plane at a given
  /// position
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the position where the measurement frame is defined
  /// @param momentum is the momentum vector (ignored)
  /// @return rotation matrix that defines the measurement frame
  RotationMatrix3D referenceFrame(const GeometryContext& gctx,
                                  const Vector3D& position,
                                  const Vector3D& momentum) const final;

  /// Return the surface type
  SurfaceType type() const override;

  /// Return method for surface normal information
  /// @note for a Cylinder a local position is always required for the normal
  /// vector
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local postion for which the normal vector is
  /// requested
  ///
  /// @return normal vector at the local position by value
  Vector3D normal(const GeometryContext& gctx,
                  const Vector2D& lposition) const final;

  /// Return method for surface normal information
  /// @note for a Cylinder a local position is always required for the normal
  /// vector
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global postion for which the normal vector is
  /// requested
  ///
  /// @return normal vector at the global position by value
  Vector3D normal(const GeometryContext& gctx,
                  const Vector3D& position) const final;

  /// Normal vector return without argument
  using Surface::normal;

  /// Return method for the rotational symmetry axis
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @return  the z-Axis of transform
  virtual Vector3D rotSymmetryAxis(const GeometryContext& gctx) const;

  /// This method returns the CylinderBounds by reference
  const CylinderBounds& bounds() const final;

  /// Local to global transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position to be transformed
  /// @param momentum is the global momentum (ignored in this operation)
  ///
  /// @return The global position by value
  Vector3D localToGlobal(const GeometryContext& gctx, const Vector2D& lposition,
                         const Vector3D& momentum) const final;

  /// Global to local transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position to be transformed
  /// @param momentum is the global momentum (ignored in this operation)
  ///
  /// @return a Result<Vector2D> which can be !ok() if the operation fails
  Result<Vector2D> globalToLocal(const GeometryContext& gctx,
                                 const Vector3D& position,
                                 const Vector3D& momentum) const final;

  /// Straight line intersection schema from position/direction
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position to start from
  /// @param direction The direction at start
  /// @param bcheck the Boundary Check
  ///
  /// If possible returns both solutions for the cylinder
  ///
  /// @return SurfaceIntersection object (contains intersection & surface)
  SurfaceIntersection intersect(const GeometryContext& gctx,
                                const Vector3D& position,
                                const Vector3D& direction,
                                const BoundaryCheck& bcheck) const final;

  /// Path correction due to incident of the track
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position as a starting point
  /// @param direction is the global momentum at the starting point
  ///
  /// @return is the correction factor due to incident
  double pathCorrection(const GeometryContext& gctx, const Vector3D& position,
                        const Vector3D& direction) const final;

  /// Return method for properly formatted output string
  std::string name() const override;

  /// Return a Polyhedron for a cylinder
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lseg Number of segments along curved lines, it represents
  /// the full 2*M_PI coverange, if lseg is set to 1 only the extrema
  /// are given
  ///
  /// @return A list of vertices and a face/facett description of it
  Polyhedron polyhedronRepresentation(const GeometryContext& gctx,
                                      size_t lseg) const override;

  /// Calculate the derivative of bound track parameters local position w.r.t.
  /// position in local 3D Cartesian coordinates
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position of the paramters in global
  ///
  /// @return Derivative of bound local position w.r.t. position in local 3D
  /// cartesian coordinates
  LocalCartesianToBoundLocalMatrix localCartesianToBoundLocalDerivative(
      const GeometryContext& gctx, const Vector3D& position) const final;

 protected:
  std::shared_ptr<const CylinderBounds> m_bounds;  //!< bounds (shared)

 private:
  /// Implementation of the intersection solver
  ///
  ///  <b>mathematical motivation:</b>
  ///
  /// The cylinder is given by :
  ///   - cylinder center: ccenter (C)
  ///   - the direction of the cylinder axis: cdirection (DZ)
  ///   - the radius r
  /// The line is given by :
  ///   - a reference position : lposition (L0)
  ///   - the line direction: ldirection (DL)
  ///   the parametric form for the line is then : L(t) = L0 + t * DL
  ///
  /// Any point P on infinite cylinder if :
  ///      ((P - C) x DZ)^2 = r^2 * DZ^2
  /// We know that DZ is a unit vector:
  ///   DZ^2 == 1
  /// When expanded with the line equation, this is  :
  ///      ((L0 - C) x DZ + t * (DL x DZ))^2 = r^2 * DZ^2
  /// which is a quadratic equation in the form (X + t * Y)^2 = d, where :
  ///  X = (L0 - C) x DZ
  ///  Y = DL x DZ
  ///  d = r^2 * (DZ)^2
  /// Now, expand the equation :
  /// t^2 * (Y . Y) + t * (2 * (X . Y)) + (X . X) - d = 0
  /// => second order equation in the form : a*t^2 + b*t + c = 0 where
  /// a = (Y . Y)
  /// b = 2 * (X . Y)
  /// c = (X . X) - d
  /// finally solve the second order equation : a*t^2 + b*t + c = 0
  /// reinsertion into the line equation.
  ///
  /// @return the quadratic equation
  detail::RealQuadraticEquation intersectionSolver(
      const Transform3D& transform, const Vector3D& position,
      const Vector3D& direction) const;
};

#include "Acts/Surfaces/detail/CylinderSurface.ipp"

}  // namespace Acts
