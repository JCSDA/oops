/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/util/ObjectCounter.h"

// --------------------------------------------------------------------------------------
// forward declarations
namespace oops {
  class GeometryData;
}

// --------------------------------------------------------------------------------------

namespace oops {

/// An explicit pseudo-diffusion smoothing operator.
class Diffusion :  private util::ObjectCounter<Diffusion>  {
 public:
  static const std::string classname() {return "oops::Diffusion";}

  /// The type of diffusion to perform (2D, 3D, etc..).
  enum class Mode {
    Joint3D,  ///< A true joint 3D diffusion (not implemented, yet)
    Split3D,  ///< A 3D diffusion where horizontal and vertical are handled separately (default)
    HorizontalOnly,  ///< Only horizontal diffusion is used
    VerticalOnly  ///< Only vertical diffusion is used
  };

  /// The numerical scheme used for vertical diffusion.
  /// Horizontal always uses the explicit scheme.
  enum class VerticalMethod {
    Explicit,  ///< Explicit iterative pseudo-diffusion (Weaver & Courtier 2001). Default.
    Implicit   ///< Implicit diffusion via per-column Cholesky/LDL^T tridiagonal solves
               ///  (Mirouze & Weaver 2010). Unconditionally stable; yields Matern-M kernels.
  };

  // derived grid geometry
  struct DerivedGeom;

  /// Initialize the geometry used by the diffusion operator.
  /// @param derivedGeom if provided, the derived geometry will not be recalculated.
  explicit Diffusion(const oops::GeometryData &,
                     const std::shared_ptr<DerivedGeom> & derivedGeom = nullptr);

  /// Calculate the derived geometry needed by Diffusion. This typically only
  /// needs to be done once for a given geometry, and can be shared among
  /// multiple Diffusion classes
  static std::shared_ptr<DerivedGeom> calculateDerivedGeom(const oops::GeometryData &);

  /// The per-node inverse control-volume area (1/area) used by the finite-volume
  /// diffusion. Exposed for verification of the area estimate on irregular meshes.
  const atlas::Field & inverseArea() const;

  /// Set the parameters used by the diffusion operator.
  /// TODO(Travis) currently this is only taking the hz and vt scales, but will
  /// be expanded to also take the precalculated diffusion constants /
  /// niter, and/or modulation fields.
  /// @param parameters contains 1 or both of the following fields:
  ///   - "hzScales" A 3D or 2D field with the horizontal length scales (units: meters).
  ///                If a 2D field is given the same scales are used on each level.
  ///   - "vtScales" A 3D field with the vertical length scales (units: number of levels)
  /// @param vtMethod one of VerticalMethod::Explicit (default) or VerticalMethod::Implicit.
  /// @param vtImplicitIterations number of implicit iterations M (must be even, >= 2).
  ///        Ignored when vtMethod is Explicit.
  void setParameters(const atlas::FieldSet & parameters,
                     VerticalMethod vtMethod = VerticalMethod::Explicit,
                     int vtImplicitIterations = 4);

  /// Perform diffusion smoothing of the input fields.
  /// If you need an operation that is self-adjoint, use `multiplySqrtAD()` and
  /// `multiplySqrtTL()` instead.
  void multiply(atlas::FieldSet &, Mode mode = Mode::Split3D) const;

  /// Perform the square root of the diffusion.
  /// (i.e. The tangent-linear, with half of the number of iterations compared with `multiply()`)
  void multiplySqrtTL(atlas::FieldSet &, Mode mode = Mode::Split3D) const;

  /// Perform adjoint of the square root of the diffusion.
  /// (i.e. The adjoint, with half the number of iterations compared with `multiply()`)
  void multiplySqrtAD(atlas::FieldSet &, Mode mode = Mode::Split3D) const;

 private:
  const oops::GeometryData & geom_;
  std::shared_ptr<DerivedGeom> derivedGeom_;

  // horizontal diffusion parameters
  int niterHz_ = -1;  // number of iterations for horizontal diffusion, or -1 if off
  int khdtLevels_;
  std::vector<std::vector<double> > khdt_;  // the horizontal diffusion constants (dim[edge][lvl] )

  // vertical diffusion parameters
  int niterVt_ = -1;  // total number of iterations M for vertical diffusion (even), or -1 if off
  atlas::Field kvdt_;  // the vertical diffusion constants (used by explicit scheme)

  // Implicit vertical diffusion state. Only populated when vtMethod_ == Implicit.
  // The tridiagonal matrix is symmetric positive-definite; we store an LDL^T
  // factorization per column, precomputed once in setParameters and reused by
  // every solve. Layout:
  //   vtImplicitD_         [ngrid, nz]    diagonal entries of D
  //   vtImplicitSubDiag_   [ngrid, nz-1]  subdiagonal of L (strictly below)
  VerticalMethod vtMethod_ = VerticalMethod::Explicit;
  atlas::Field vtImplicitD_;
  atlas::Field vtImplicitSubDiag_;

  // Square-root building blocks. Each routine applies M/2 iterations (the
  // "square root" of the full diffusion); multiply() composes two of them
  // around the horizontal pass to deliver the full M-iteration operator.
  void multiplyHzTL(atlas::Field &) const;
  void multiplyHzAD(atlas::Field &) const;
  void multiplyVtExplicitTL(atlas::Field &) const;
  void multiplyVtExplicitAD(atlas::Field &) const;
  // Implicit vertical solve is exactly self-adjoint (the discrete tridiagonal
  // matrix is symmetric by construction), so a single routine serves as both
  // square-root TL and AD.
  void multiplyVtImplicit(atlas::Field &) const;
};

// --------------------------------------------------------------------------------------

}  // namespace oops
