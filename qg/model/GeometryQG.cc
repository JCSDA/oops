/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "model/GeometryQG.h"
#include "model/QgFortran.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
GeometryQG::GeometryQG(const GeometryQgParameters & params,
                       const eckit::mpi::Comm & comm) : comm_(comm), levs_(0) {
  ASSERT(comm_.size() == 1);
  qg_geom_setup_f90(keyGeom_, params.toConfiguration());

  int nx = 0;
  int ny = 0;
  int nz;
  double deltax;
  double deltay;
  qg_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay);
  levs_ = nz;

  // Set ATLAS lon/lat field
  atlasFieldSet_.reset(new atlas::FieldSet());
  qg_geom_set_atlas_lonlat_f90(keyGeom_, atlasFieldSet_->get());
  atlas::Field atlasField = atlasFieldSet_->field("lonlat");

  // Create ATLAS function space
  atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(atlasField));

  // Set ATLAS function space pointer in Fortran
  qg_geom_set_atlas_functionspace_pointer_f90(keyGeom_, atlasFunctionSpace_->get());

  // Fill ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  qg_geom_fill_atlas_fieldset_f90(keyGeom_, atlasFieldSet_->get());
}
// -----------------------------------------------------------------------------
GeometryQG::GeometryQG(const GeometryQG & other) : comm_(other.comm_), levs_(other.levs_) {
  ASSERT(comm_.size() == 1);
  qg_geom_clone_f90(keyGeom_, other.keyGeom_);

  // Copy ATLAS function space
  atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(
                            other.atlasFunctionSpace_->lonlat()));

  // Set ATLAS function space pointer in Fortran
  qg_geom_set_atlas_functionspace_pointer_f90(keyGeom_, atlasFunctionSpace_.get()->get());

  // Copy ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  for (int jfield = 0; jfield < other.atlasFieldSet_->size(); ++jfield) {
    atlas::Field atlasField = other.atlasFieldSet_->field(jfield);
    atlasFieldSet_->add(atlasField);
  }
}
// -----------------------------------------------------------------------------
GeometryQG::~GeometryQG() {
  qg_geom_delete_f90(keyGeom_);
}
// -----------------------------------------------------------------------------
GeometryQGIterator GeometryQG::begin() const {
  return GeometryQGIterator(*this);
}
// -----------------------------------------------------------------------------
GeometryQGIterator GeometryQG::end() const {
  int nx = 0;
  int ny = 0;
  int nz = 1;
  double deltax;
  double deltay;
  qg_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay);
  return GeometryQGIterator(*this, nx*ny+1);
}
// -------------------------------------------------------------------------------------------------
void GeometryQG::latlon(std::vector<double> & lats, std::vector<double> & lons, const bool) const {
  const auto lonlat = atlas::array::make_view<double, 2>(atlasFunctionSpace_->lonlat());
  const size_t npts = atlasFunctionSpace_->size();
  lats.resize(npts);
  lons.resize(npts);
  for (size_t jj = 0; jj < npts; ++jj) {
    lats[jj] = lonlat(jj, 1);
    lons[jj] = lonlat(jj, 0);
  }
}
// -------------------------------------------------------------------------------------------------
std::vector<double> GeometryQG::verticalCoord(std::string & vcUnits) const {
  // returns vertical coordinate in untis of vcUnits
  int nx = 0;
  int ny = 0;
  int nz = 1;
  double deltax;
  double deltay;
  qg_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay);
  std::vector<double> vc(nz);
  if (vcUnits == "levels") {
    for (int i=0; i < nz; ++i) {vc[i]=i+1;}
  } else {
    std::stringstream errorMsg;
    errorMsg << "Uknown vertical coordinate unit " << vcUnits << std::endl;
    ABORT(errorMsg.str());
  }
  oops::Log::debug() << "QG vert coord: " << vc << std::endl;
  return vc;
}
// -------------------------------------------------------------------------------------------------
std::vector<size_t> GeometryQG::variableSizes(const oops::Variables & vars) const {
  std::vector<size_t> sizes(vars.size(), levs_);
  return sizes;
}
// -----------------------------------------------------------------------------
void GeometryQG::print(std::ostream & os) const {
  int nx;
  int ny;
  int nz;
  double deltax;
  double deltay;
  qg_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay);
  os << "Geometry:" << std::endl;
  os << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
  os << "deltax = " << deltax << ", deltay = " << deltay;
}
// -----------------------------------------------------------------------------
}  // namespace qg
