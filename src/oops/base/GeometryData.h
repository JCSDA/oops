/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"

namespace eckit {
namespace mpi {
class Comm;
}  // namespace mpi
}  // namespace eckit

namespace stripack {
class Triangulation;
}  // namespace stripack

namespace oops {

// -----------------------------------------------------------------------------

class GeometryData {
 public:
  GeometryData(const atlas::FunctionSpace &, const atlas::FieldSet &,
               const bool, const eckit::mpi::Comm &);

  ~GeometryData();  // defined in .cc file, where stripack::Triangulation is complete

// Local tree requires lats and lons with halo
  void setLocalTree(const std::vector<double> &, const std::vector<double> &);
// Global tree requires lats and lons without halo
  void setGlobalTree(const std::vector<double> &, const std::vector<double> &);

  GeometryData(const GeometryData &) = delete;
  GeometryData & operator=(const GeometryData &) = delete;

  int closestTask(const double, const double) const;
  atlas::util::KDTree<size_t>::ValueList closestPoints(const double, const double, const int) const;

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
  atlas::FunctionSpace fspace_;
  atlas::FieldSet fset_;
  const eckit::mpi::Comm * comm_;
  bool topdown_;

  const atlas::Geometry earth_;
  atlas::util::IndexKDTree localTree_;
  atlas::util::IndexKDTree globalTree_;
  bool loctree_;
  bool glotree_;

  const atlas::Geometry unitsphere_;
  std::vector<double> lats_;
  std::vector<double> lons_;
  // Triangulation is a bit expensive (and not valid for models like L95), so compute on demand
  mutable std::unique_ptr<const stripack::Triangulation> triangulation_;
};

// -----------------------------------------------------------------------------

}  // namespace oops
