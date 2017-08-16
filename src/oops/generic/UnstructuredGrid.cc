/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */
#include "eckit/config/Configuration.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/generic/unstructured_grid_f.h"
#include "oops/generic/UnstructuredGrid.h"
namespace oops {

// -----------------------------------------------------------------------------
UnstructuredGrid::UnstructuredGrid() : keyUGrid_(0) {
  create_ug_f90(keyUGrid_);
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
void UnstructuredGrid::dirac(const eckit::Configuration & conf) {
  const eckit::Configuration * fconf = &conf;
  setup_dirac_f90(keyUGrid_, &fconf);
}
// -----------------------------------------------------------------------------
double UnstructuredGrid::dot_product_with(const UnstructuredGrid & other) const {
  double zz = 0.0;
//  const int keyOvecOther = other.keyUGrid_;
//  dotprod_f90(keyUGrid_, keyOvecOther, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void UnstructuredGrid::print(std::ostream & os) const {
  os << " UnstructuredGrid: print not implemented yet.";
}
// -----------------------------------------------------------------------------

}  // namespace oops
