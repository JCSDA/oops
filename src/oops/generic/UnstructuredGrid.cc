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
std::vector<double> UnstructuredGrid::getLats() {
  int ncols;
  get_ncols_f90(keyUGrid_, ncols);
  std::vector<double> lats(ncols);
  get_lats_f90(keyUGrid_, ncols, &lats[0]);
  return lats;
}
// -----------------------------------------------------------------------------
std::vector<double> UnstructuredGrid::getLons() {
  int ncols;
  get_ncols_f90(keyUGrid_, ncols);
  std::vector<double> lons(ncols);
  get_lons_f90(keyUGrid_, ncols, &lons[0]);
  return lons;
}
// -----------------------------------------------------------------------------
std::vector<double> UnstructuredGrid::getAreas() {
  int ncols;
  get_ncols_f90(keyUGrid_, ncols);
  std::vector<double> areas(ncols);
  get_areas_f90(keyUGrid_, ncols, &areas[0]);
  return areas;
}
// -----------------------------------------------------------------------------
std::vector<double> UnstructuredGrid::getVunit() {
  int nlevs;
  get_nlevs_f90(keyUGrid_, nlevs);
  std::vector<double> vunit(nlevs);
  get_vunit_f90(keyUGrid_, nlevs, &vunit[0]);
  return vunit;
}
// -----------------------------------------------------------------------------
std::vector<int> UnstructuredGrid::getMask(const int & ilev) {
  int ncols;
  get_ncols_f90(keyUGrid_, ncols);
  std::vector<int> mask(ncols);
  get_mask_f90(keyUGrid_, ncols, ilev, &mask[0]);
  return mask;
}
// -----------------------------------------------------------------------------
std::vector<int> UnstructuredGrid::getGlbInd() {
  int ncols;
  get_ncols_f90(keyUGrid_, ncols);
  std::vector<int> glbind(ncols);
  get_glbind_f90(keyUGrid_, ncols, &glbind[0]);
  return glbind;
}
// -----------------------------------------------------------------------------
int UnstructuredGrid::getNvar() {
  int nvar;
  get_nvar_f90(keyUGrid_, nvar);
  return nvar;
}
// -----------------------------------------------------------------------------
std::vector<double> UnstructuredGrid::getData() {
  int nlevs;
  get_nlevs_f90(keyUGrid_, nlevs);
  int ncols;
  get_ncols_f90(keyUGrid_, ncols);
  int nvar;
  get_nvar_f90(keyUGrid_, nvar);
  std::vector<double> data(nlevs*ncols*nvar);
  get_data_f90(keyUGrid_, nlevs*ncols*nvar, &data[0]);
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
