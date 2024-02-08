/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/mesh.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"

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
               const bool, const eckit::mpi::Comm &);

  ~GeometryData() = default;

// Global tree requires lats and lons without halo
  void setGlobalTree(const std::vector<double> &, const std::vector<double> &);

  GeometryData(const GeometryData &) = delete;
  GeometryData & operator=(const GeometryData &) = delete;

  int closestTask(const double, const double) const;

  /// Identifies the three model grid points defining the triangle containing (lat,lon).
  ///
  /// Returns true if such a triangle is found; false if not.
  bool containingTriangleAndBarycentricCoords(double lat, double lon,
      std::array<int, 3> & indices, std::array<double, 3> & barycentricCoords) const;

// Accessors
  const atlas::FunctionSpace & functionSpace() const {return fspace_;}
  const atlas::FieldSet & fieldSet() const {return fset_;}
  const eckit::mpi::Comm & comm() const {return *comm_;}

  bool has(const std::string & name) const {return fset_.has(name);}
  const atlas::Field & getField(const std::string & name) const {return fset_.field(name);}
  bool levelsAreTopDown() const {return topdown_;}

 private:
  void setMeshAndTriangulation();
  void setLocalTree();

  atlas::FunctionSpace fspace_;
  atlas::FieldSet fset_;
  const eckit::mpi::Comm * comm_;
  bool topdown_;

  atlas::Mesh mesh_;
  std::vector<bool> firstTriangulationOfQuadsIsDelaunay_;
  const atlas::Geometry earth_;
  atlas::util::IndexKDTree globalNodeTree_;  // JEDI grid nodes = model cell-centers
  atlas::util::IndexKDTree localCellCenterTree_;  // JEDI cell centers
};

// -----------------------------------------------------------------------------

}  // namespace oops
