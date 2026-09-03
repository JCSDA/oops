/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <array>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/mesh.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"

#include "oops/util/FunctionSpaceHelpers.h"

namespace eckit {
namespace mpi {
class Comm;
}  // namespace mpi
}  // namespace eckit

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Mesh-triangulation data structures used by the oops::UnstructuredInterpolator
///        for locating target points within the grid, derived from a source FunctionSpace.
///
/// Construction performs MPI-collective work (mesh generation, halo build, gather /
/// broadcast for custom partitions), so it must be reached uniformly across all ranks
/// of \p comm.
class MeshTriangulation {
 public:
  MeshTriangulation(const atlas::FunctionSpace &, const eckit::mpi::Comm &);

  ~MeshTriangulation() = default;

  MeshTriangulation(const MeshTriangulation &) = delete;
  MeshTriangulation & operator=(const MeshTriangulation &) = delete;

  /// Identifies the three model grid points defining the triangle containing (lat,lon).
  ///
  /// Returns true if such a triangle is found; false if not.
  bool containingTriangleAndBarycentricCoords(double lat, double lon,
      std::array<int, 3> & indices, std::array<double, 3> & barycentricCoords) const;

 private:
  void setMeshAndTriangulation(const atlas::FunctionSpace &, const eckit::mpi::Comm &);
  void setLocalTree();

  atlas::Mesh mesh_;
  std::vector<bool> firstTriangulationOfQuadsIsDelaunay_;
  const atlas::Geometry earth_;
  atlas::util::IndexKDTree localCellCenterTree_;  // JEDI cell centers

  // When the FunctionSpace is StructuredColumns, we'll need to map the Mesh indices used in the
  // interpolation stencil computation onto the FunctionSpace indices used to read FieldSet data.
  bool is_atlas_structured_columns_ = false;
  util::StructuredMeshToStructuredColumnsIndexMap indexMapper_;
  // When the grid is regular (in either StructuredColumns or NodeColumns), we need to look through
  // more triangles at the poles to account for grid deformation.
  int regular_grid_nx_ = -1;
};

// -----------------------------------------------------------------------------

/// \brief Returns the MeshTriangulation for the source geometry described by
///        (\p fspace, \p comm), building and caching it on first request.
///
/// The cache is held in a process-local registry keyed by grid identity (communicator
/// name + grid UID + lon/lat hash), so repeated interpolators on the same source
/// geometry share a single mesh/tree build. Because construction is MPI-collective,
/// this must be called uniformly across all ranks of \p comm.
const MeshTriangulation & getMeshTriangulation(
    const atlas::FunctionSpace & fspace, const eckit::mpi::Comm & comm);

// -----------------------------------------------------------------------------

}  // namespace oops
