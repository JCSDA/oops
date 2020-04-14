/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "eckit/config/Configuration.h"

#include "model/GeometryQG.h"
#include "model/QgFortran.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
GeometryQG::GeometryQG(const eckit::Configuration & conf,
                       const eckit::mpi::Comm & comm) : comm_(comm) {
  qg_geom_setup_f90(keyGeom_, conf);

  // Create ATLAS grid configuration
  const atlas::util::Config atlasConfig;
  qg_geom_create_atlas_grid_conf_f90(keyGeom_, atlasConfig);

  // Create ATLAS grid
  atlas::StructuredGrid atlasStructuredGrid(atlasConfig);

  // Create ATLAS function space
  atlasFunctionSpace_.reset(new atlas::functionspace::StructuredColumns(atlasStructuredGrid,
                            atlas::option::halo(0)));

  // Set ATLAS function space pointer in Fortran
  qg_geom_set_atlas_functionspace_pointer_f90(keyGeom_, atlasFunctionSpace_.get()->get());

  // Fill ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  qg_geom_fill_atlas_fieldset_f90(keyGeom_, atlasFieldSet_.get()->get());

  // Set ATLAS fieldset pointer in Fortran
  qg_geom_set_atlas_fieldset_pointer_f90(keyGeom_, atlasFieldSet_.get()->get());
}
// -----------------------------------------------------------------------------
GeometryQG::GeometryQG(const GeometryQG & other) : comm_(other.comm_) {
  qg_geom_clone_f90(keyGeom_, other.keyGeom_);

  // Copy ATLAS function space
  atlasFunctionSpace_.reset(new atlas::functionspace::StructuredColumns(
                            other.atlasFunctionSpace_->grid(), atlas::option::halo(0)));

  // Set ATLAS function space pointer in Fortran
  qg_geom_set_atlas_functionspace_pointer_f90(keyGeom_, atlasFunctionSpace_.get()->get());

  // Copy ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  for (int jfield = 0; jfield < other.atlasFieldSet_->size(); ++jfield) {
    atlas::Field atlasField = other.atlasFieldSet_->field(jfield);
    atlasFieldSet_->add(atlasField);
  }

  // Set ATLAS fieldset pointer in Fortran
  qg_geom_set_atlas_fieldset_pointer_f90(keyGeom_, atlasFieldSet_.get()->get());
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
  int nz;
  double deltax;
  double deltay;
  qg_geom_info_f90(keyGeom_, nx, ny, nz, deltax, deltay);
  return GeometryQGIterator(*this, nx*ny+1);
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
  os << "deltax = " << deltax << ", deltay = " << deltay << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace qg
