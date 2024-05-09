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
#include "atlas/mesh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
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

  // Query coordinates to build FunctionSpace
  atlas::FieldSet fieldSet;
  qg_geom_set_lonlat_f90(keyGeom_, fieldSet.get());
  const auto lonlats = atlas::array::make_view<double, 2>(fieldSet.field("lonlat"));

  // QG grid is similar to an atlas longitude-shifted RegularLonLatGrid, or "Slon<NLON>x<NLAT>",
  // except that the mapping of points in the latitude direction isn't obviously reproducible.
  // So, import QG grid into atlas via MeshBuilder:
  const int npoints = nx * ny;
  const int ncells = nx * (ny - 1);
  std::vector<double> lons(npoints);
  std::vector<double> lats(npoints);
  std::vector<int> ghosts(npoints, 0);
  std::vector<atlas::gidx_t> global_indices(npoints);
  std::vector<atlas::idx_t> remote_indices(npoints);
  atlas::idx_t remote_index_base = 1;
  std::vector<int> partitions(npoints, 0);
  std::vector<std::array<atlas::gidx_t, 3>> tri_boundary_nodes{};
  std::vector<atlas::gidx_t> tri_global_indices{};
  std::vector<std::array<atlas::gidx_t, 4>> quad_boundary_nodes(ncells);
  std::vector<atlas::gidx_t> quad_global_indices(ncells);

  int jpoint = 0;
  int jcell = 0;
  for (int jy = 0; jy < ny; ++jy) {
    for (int jx = 0; jx < nx; ++jx) {
      lons[jpoint] = lonlats(jpoint, 0);
      lats[jpoint] = lonlats(jpoint, 1);
      const int jindex = jpoint + 1;
      global_indices[jpoint] = jindex;
      remote_indices[jpoint] = jindex;
      if (jy != ny-1) {
        if (jx != nx-1) {
          quad_boundary_nodes[jcell] = {jindex, jindex + 1, jindex + nx + 1, jindex + nx};
        } else {
          quad_boundary_nodes[jcell] = {jindex, jindex + 1 - nx, jindex + 1, jindex + nx};
        }
        quad_global_indices[jcell] = jcell + 1;
        jcell++;
      }
      jpoint++;
    }
  }
  ASSERT(jpoint == npoints);
  ASSERT(jcell == ncells);

  eckit::LocalConfiguration config{};
  config.set("mpi_comm", comm_.name());
  const atlas::mesh::MeshBuilder mesh_builder{};
  const atlas::Mesh mesh = mesh_builder(lons, lats, ghosts, global_indices,
                                        remote_indices, remote_index_base, partitions,
                                        tri_boundary_nodes, tri_global_indices,
                                        quad_boundary_nodes, quad_global_indices, config);
  // Don't call build_halo: QG is a serial model so there are no halos
  functionSpace_ = atlas::functionspace::NodeColumns(mesh, config);

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
