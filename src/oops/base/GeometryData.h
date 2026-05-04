/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
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

class GeometryData {
 public:
  GeometryData(const atlas::FunctionSpace &, const atlas::FieldSet &,
               const bool, const eckit::mpi::Comm &,
               const bool & initMeshAndTrees = true);

  ~GeometryData() = default;

  GeometryData(const GeometryData &) = delete;
  GeometryData & operator=(const GeometryData &) = delete;

  int closestTask(const double, const double) const;

  /// Returns the task-local index into the FunctionSpace of the globally-nearest
  /// grid point to (lat, lon). Aborts if the tree is empty or if the nearest point
  /// is not owned by this MPI rank. It is assumed that the closestTask method will
  /// have been previously used, so that the call to this method is placed on the
  /// correct MPI task.
  /// Note that `radius` is a chord distance in meters.
  std::optional<int> closestPointWithinRadius(double lat, double lon, double radius) const;

  /// Identifies the three model grid points defining the triangle containing (lat,lon).
  ///
  /// Returns true if such a triangle is found; false if not.
  bool containingTriangleAndBarycentricCoords(double lat, double lon,
      std::array<int, 3> & indices, std::array<double, 3> & barycentricCoords) const;

  const atlas::FunctionSpace & functionSpace() const {return fspace_;}
  const atlas::FieldSet & fieldSet() const {return fset_;}
  const eckit::mpi::Comm & comm() const {return comm_;}

  bool has(const std::string & name) const {return fset_.has(name);}
  const atlas::Field & getField(const std::string & name) const {return fset_.field(name);}
  bool levelsAreTopDown() const {return topdown_;}

 private:
  void setGlobalTree();
  void setMeshAndTriangulation();
  void setLocalTree();

  atlas::FunctionSpace fspace_;
  atlas::FieldSet fset_;
  const eckit::mpi::Comm & comm_;
  bool topdown_;

  atlas::Mesh mesh_;
  std::vector<bool> firstTriangulationOfQuadsIsDelaunay_;
  const atlas::Geometry earth_;
  // JEDI grid nodes = model cell-centers
  atlas::util::KDTree<std::pair<atlas::idx_t, atlas::idx_t>> globalNodeTree_;
  // JEDI cell centers
  atlas::util::IndexKDTree localCellCenterTree_;

  // When the FunctionSpace is StructuredColumns, we'll need to map the Mesh indices used in the
  // interpolation stencil computation onto the FunctionSpace indices used to read FieldSet data.
  bool is_atlas_structured_columns_ = false;
  util::StructuredMeshToStructuredColumnsIndexMap indexMapper_;
  // When the grid is regular (in either StructuredColumns or NodeColumns), we need to look through
  // more triangles at the poles to account for grid deformation.
  int regular_grid_nx_ = -1;
};

// -----------------------------------------------------------------------------

}  // namespace oops
