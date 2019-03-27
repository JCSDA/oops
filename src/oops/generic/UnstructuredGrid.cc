/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "eckit/config/Configuration.h"
#include "oops/generic/unstructured_grid_f.h"
#include "oops/generic/UnstructuredGrid.h"

namespace oops {

// -----------------------------------------------------------------------------
UnstructuredGrid::UnstructuredGrid(const int & colocated, const int & nts) : keyUGrid_(0),
  colocated_(colocated), nts_(nts) {
  create_ug_f90(keyUGrid_, colocated_, nts_);
}
// -----------------------------------------------------------------------------
UnstructuredGrid::~UnstructuredGrid() {
  delete_ug_f90(keyUGrid_);
}
// -----------------------------------------------------------------------------
void UnstructuredGrid::zero() {
}
// -----------------------------------------------------------------------------
void UnstructuredGrid::random() {
}
// -----------------------------------------------------------------------------
double UnstructuredGrid::dot_product_with(const UnstructuredGrid & other) const {
  double zz = 0.0;
  return zz;
}
// -----------------------------------------------------------------------------
void UnstructuredGrid::print(std::ostream & os) const {
  os << " UnstructuredGrid: print not implemented yet.";
}
// -----------------------------------------------------------------------------

}  // namespace oops
