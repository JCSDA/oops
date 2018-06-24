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
UnstructuredGrid::UnstructuredGrid() : keyUGrid_(0) {
  create_ug_f90(keyUGrid_);
}
// -----------------------------------------------------------------------------
UnstructuredGrid::~UnstructuredGrid() {
  delete_ug_f90(keyUGrid_);
}
// -----------------------------------------------------------------------------
int UnstructuredGrid::getSize(const int & ind) {
  int isize;
  get_size_f90(keyUGrid_, ind, isize);
  return isize;
}
// -----------------------------------------------------------------------------
std::vector<double> UnstructuredGrid::getLon() {
  int nmga;
  get_size_f90(keyUGrid_, 1, nmga);
  std::vector<double> lon(nmga);
  get_lon_f90(keyUGrid_, nmga, &lon[0]);
  return lon;
}
// -----------------------------------------------------------------------------
std::vector<double> UnstructuredGrid::getLat() {
  int nmga;
  get_size_f90(keyUGrid_, 1, nmga);
  std::vector<double> lat(nmga);
  get_lat_f90(keyUGrid_, nmga, &lat[0]);
  return lat;
}
// -----------------------------------------------------------------------------
std::vector<double> UnstructuredGrid::getArea() {
  int nmga;
  get_size_f90(keyUGrid_, 1, nmga);
  std::vector<double> area(nmga);
  get_area_f90(keyUGrid_, nmga, &area[0]);
  return area;
}
// -----------------------------------------------------------------------------
std::vector<double> UnstructuredGrid::getVunit() {
  int nmga;
  get_size_f90(keyUGrid_, 1, nmga);
  int nl0;
  get_size_f90(keyUGrid_, 2, nl0);
  std::vector<double> vunit(nmga*nl0);
  get_vunit_f90(keyUGrid_, nmga*nl0, &vunit[0]);
  return vunit;
}
// -----------------------------------------------------------------------------
std::vector<int> UnstructuredGrid::getImask() {
  int nmga;
  get_size_f90(keyUGrid_, 1, nmga);
  int nl0;
  get_size_f90(keyUGrid_, 2, nl0);
  std::vector<int> imask(nmga*nl0);
  get_imask_f90(keyUGrid_, nmga*nl0, &imask[0]);
  return imask;
}
// -----------------------------------------------------------------------------
std::vector<double> UnstructuredGrid::getData() {
  int nmga;
  get_size_f90(keyUGrid_, 1, nmga);
  int nl0;
  get_size_f90(keyUGrid_, 2, nl0);
  int nv;
  get_size_f90(keyUGrid_, 3, nv);
  int nts;
  get_size_f90(keyUGrid_, 4, nts);
  std::vector<double> data(nmga*nl0*nv*nts);
  get_data_f90(keyUGrid_, nmga*nl0*nv*nts, &data[0]);
  return data;
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
