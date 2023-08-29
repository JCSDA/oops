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
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "model/GeometryQG.h"
#include "model/QgFortran.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
GeometryQG::GeometryQG(const eckit::Configuration & conf,
                       const eckit::mpi::Comm & comm) : comm_(comm), levs_(0) {
  ASSERT(comm_.size() == 1);
  qg_geom_setup_f90(keyGeom_, conf);

  int nx = 0;
  int ny = 0;
  int nz;
  double deltax;
  double deltay;
  qg_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay);
  levs_ = nz;

  // Create function space
  atlas::FieldSet fieldSet;
  qg_geom_set_lonlat_f90(keyGeom_, fieldSet.get());
  const atlas::Field & field = fieldSet.field("lonlat");
  functionSpace_ = atlas::functionspace::PointCloud(field);

  // Set function space pointer in Fortran
  qg_geom_set_functionspace_pointer_f90(keyGeom_, functionSpace_.get());

  // Fill geometry fields
  fields_ = atlas::FieldSet();
  qg_geom_fill_geometry_fields_f90(keyGeom_, fields_.get());
}
// -----------------------------------------------------------------------------
GeometryQG::GeometryQG(const GeometryQG & other) : comm_(other.comm_), levs_(other.levs_) {
  ASSERT(comm_.size() == 1);
  qg_geom_clone_f90(keyGeom_, other.keyGeom_);

  // Copy function space
  functionSpace_ = other.functionSpace_;

  // Set function space pointer in Fortran
  qg_geom_set_functionspace_pointer_f90(keyGeom_, functionSpace_.get());

  // Copy geometry fields
  fields_ = atlas::FieldSet();
  for (auto & field : other.fields_) {
    fields_.add(field);
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
  const auto lonlat = atlas::array::make_view<double, 2>(functionSpace_.lonlat());
  const size_t npts = functionSpace_.size();
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
