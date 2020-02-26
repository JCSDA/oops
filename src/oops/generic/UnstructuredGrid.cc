/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <iomanip>
#include <sstream>
#include <utility>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/util/Config.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Variables.h"
#include "oops/generic/unstructured_grid_f.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/util/DateTime.h"

namespace oops {

// -----------------------------------------------------------------------------
UnstructuredGrid::UnstructuredGrid(const int & colocated, const int & nts) :
  keyUGrid_(0), atlasFunctionSpace_(), atlasFieldSet_() {
  create_ug_f90(keyUGrid_, colocated, nts);
}
// -----------------------------------------------------------------------------
UnstructuredGrid::UnstructuredGrid(UnstructuredGrid & other) :
  keyUGrid_(other.keyUGrid_), atlasFunctionSpace_(std::move(other.atlasFunctionSpace_)),
  atlasFieldSet_(std::move(other.atlasFieldSet_)) {
  other.keyUGrid_ = 0;
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
void UnstructuredGrid::defineGeometry() {
  // Create ATLAS grid configuration
  const atlas::util::Config atlasConfig;
  const eckit::Configuration * fconf = &atlasConfig;
  ug_create_atlas_grid_conf_f90(keyUGrid_, &fconf);
  atlas::UnstructuredGrid atlasUnstructuredGrid(atlasConfig);

  // Create mesh
  atlas::MeshGenerator atlasMeshGenerator("no_connectivity");
  atlas::Mesh atlasMesh = atlasMeshGenerator.generate(atlasUnstructuredGrid);

  // Create ATLAS function space
  atlasFunctionSpace_.reset(new atlas::functionspace::NodeColumns(atlasMesh,
                            atlas::option::halo(0)));

  // Set ATLAS function space pointer in Fortran
  ug_set_atlas_functionspace_pointer_f90(keyUGrid_, atlasFunctionSpace_->get());

  // Fill ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  ug_fill_atlas_fieldset_f90(keyUGrid_, atlasFieldSet_->get());

  // Set ATLAS fieldset in Fortran
  ug_set_atlas_fieldset_pointer_f90(keyUGrid_, atlasFieldSet_->get());
}
// -----------------------------------------------------------------------------
void UnstructuredGrid::defineGrids(std::vector<eckit::LocalConfiguration> & grids) const {
  int ngrid;
  ug_get_ngrid_f90(keyUGrid_, ngrid);
  for (int jgrid = 0; jgrid < ngrid; ++jgrid) {
    int nl, nv, nts;
    ug_get_dims_f90(keyUGrid_, jgrid, nl, nv, nts);
    std::vector<std::string> variables;
    for (int iv = 0; iv < nv; ++iv) {
      std::ostringstream ss;
      ss << std::setw(2) << std::setfill('0') << iv+1;
      variables.push_back("var_" + ss.str());
    }
    std::vector<std::string> timeslots;
    for (int its = 0; its < nts; ++its) {
      std::ostringstream ss;
      ss << std::setw(2) << std::setfill('0') << its+1;
      timeslots.push_back(ss.str());
    }
    eckit::LocalConfiguration grid;
    grid.set("grid_index", jgrid);
    grid.set("nl", nl);
    grid.set("variables", variables);
    grid.set("timeslots", timeslots);
    grid.set("lev2d", "first");
    grids.push_back(grid);
  }
}
// -----------------------------------------------------------------------------
void UnstructuredGrid::setAtlas(atlas::FieldSet * atlasFieldSet) const {
  ug_set_atlas_f90(keyUGrid_, atlasFieldSet->get());
}
// -----------------------------------------------------------------------------
void UnstructuredGrid::toAtlas(atlas::FieldSet * atlasFieldSet) const {
  ug_to_atlas_f90(keyUGrid_, atlasFieldSet->get());
}
// -----------------------------------------------------------------------------
void UnstructuredGrid::fromAtlas(atlas::FieldSet * atlasFieldSet) {
  ug_from_atlas_f90(keyUGrid_, atlasFieldSet->get());
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
