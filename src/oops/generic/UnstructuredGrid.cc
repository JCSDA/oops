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
std::vector<double> UnstructuredGrid::getLevs() {
  int nlevs;
  get_nlevs_f90(keyUGrid_, nlevs);
  std::vector<double> levs(nlevs);
  get_levs_f90(keyUGrid_, nlevs, &levs[0]);
  return levs;
}
// -----------------------------------------------------------------------------
std::vector<int> UnstructuredGrid::getMask3d(const int & ilev) {
  int ncols;
  get_ncols_f90(keyUGrid_, ncols);
  std::vector<int> mask3d(ncols);
  get_mask3d_f90(keyUGrid_, ncols, ilev, &mask3d[0]);
  return mask3d;
}
// -----------------------------------------------------------------------------
std::vector<int> UnstructuredGrid::getMask2d() {
  int ncols;
  get_ncols_f90(keyUGrid_, ncols);
  std::vector<int> mask2d(ncols);
  get_mask2d_f90(keyUGrid_, ncols, &mask2d[0]);
  return mask2d;
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
void UnstructuredGrid::zero() {
}
// -----------------------------------------------------------------------------
void UnstructuredGrid::random() {
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
