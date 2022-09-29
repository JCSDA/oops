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
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"

#include "oops/mpi/mpi.h"

namespace oops {

// -----------------------------------------------------------------------------

class GeometryData {
 public:
  GeometryData(const atlas::FunctionSpace &, const atlas::FieldSet &,
               const bool, const eckit::mpi::Comm &);

// Local tree requires lats and lons with halo
  void setLocalTree(const std::vector<double> &, const std::vector<double> &);
// Global tree requires lats and lons without halo
  void setGlobalTree(const std::vector<double> &, const std::vector<double> &);

  GeometryData(const GeometryData &) = delete;
  GeometryData & operator=(const GeometryData &) = delete;

  int closestTask(const double, const double) const;
  atlas::util::KDTree<size_t>::ValueList closestPoints(const double, const double, const int) const;

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
};

// -----------------------------------------------------------------------------

}  // namespace oops
