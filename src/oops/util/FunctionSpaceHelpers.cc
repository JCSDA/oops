/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/FunctionSpaceHelpers.h"

#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/abor1_cpp.h"

namespace util {

// -----------------------------------------------------------------------------

void setupFunctionSpace(const eckit::mpi::Comm & comm,
    const eckit::Configuration & config,
    atlas::Grid & grid,
    atlas::grid::Partitioner & partitioner,
    atlas::Mesh & mesh,
    atlas::FunctionSpace & functionSpace,
    atlas::FieldSet & fieldset) {
  // Set up atlas MPI
  eckit::mpi::setCommDefault(comm.name().c_str());

  // Set up Grid
  const eckit::LocalConfiguration confGrid(config, "grid");
  grid = atlas::Grid(confGrid);
  // grid.type() doesn't report if the grid is regional, so we look at the config directly
  const std::string & gridType = confGrid.getString("type", "no_type");

  // Set up Partitioner
  const bool noPointOnLastTask = config.getBool("no point on last task", false);
  const std::string partitionerName = config.getString("partitioner", "equal_regions");
  if (noPointOnLastTask && (comm.size() > 1)) {
    partitioner = atlas::grid::Partitioner(partitionerName, comm.size()-1);
  } else {
    partitioner = atlas::grid::Partitioner(partitionerName);
  }

  // Set up FunctionSpace (and Mesh)
  const std::string functionSpaceName = config.getString("function space");
  if (functionSpaceName == "StructuredColumns") {
    ASSERT(gridType != "unstructured");

    const size_t halo = config.getUnsigned("halo", 0);
    functionSpace = atlas::functionspace::StructuredColumns(grid, partitioner,
                                                            atlas::option::halo(halo));

    // Bugfix for regional grids
    // It seems that the content of lonlat for a regional function space is actually the xy
    // coordinates. The routine to compute distances on the sphere was complaining about
    // impossible lon/lat values...
    if (gridType == "regional") {
      auto lonlat = atlas::array::make_view<double, 2>(functionSpace.lonlat());
      double lonlatPoint[] = {0, 0};
      const atlas::functionspace::StructuredColumns fs(functionSpace);
      const atlas::StructuredGrid grid = fs.grid();
      const auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
      const auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
      for (int jj = 0; jj < fs.size(); ++jj) {
        grid.lonlat(view_i(jj)-1, view_j(jj)-1, lonlatPoint);
        lonlat(jj, 0) = lonlatPoint[0];
        lonlat(jj, 1) = lonlatPoint[1];
      }
    }

    mesh = atlas::MeshGenerator("structured").generate(grid, partitioner);
  } else if (functionSpaceName == "NodeColumns") {
    if (grid.name().compare(0, 2, std::string{"CS"}) == 0) {
      // NodeColumns from a CubedSphere grid/mesh
      mesh = atlas::MeshGenerator("cubedsphere_dual").generate(grid, partitioner);
      functionSpace = atlas::functionspace::CubedSphereNodeColumns(mesh);
    } else {
      // NodeColumns from an unstructured grid after triangulation
      // Regular or Structured grids could be supported with extra code
      mesh = atlas::MeshGenerator("delaunay").generate(grid, partitioner);
      functionSpace = atlas::functionspace::NodeColumns(mesh);
    }
  } else if (functionSpaceName == "PointCloud") {
    ABORT(functionSpaceName + " function space not supported");
  } else {
    ABORT(functionSpaceName + " function space not supported yet");
  }

  // Owned points mask, typically the inverse of atlas' ghost mask
  atlas::Field owned = functionSpace.createField<int>(atlas::option::name("owned")
                                                      | atlas::option::levels(1));
  auto ownedView = atlas::array::make_view<int, 2>(owned);
  auto ghostView = atlas::array::make_view<int, 1>(functionSpace.ghost());
  for (atlas::idx_t i = 0; i < ghostView.shape(0); ++i) {
    ownedView(i, 0) = (ghostView(i) > 0 ? 0 : 1);
  }

  // Owned points -- special case for lat/lon grids
  // TODO(algo): this if-block can likely be generalized to other grid types (Gaussian) where
  //             atlas also has "self halos" at the grid seams, perhaps by using an algorithm
  //             to identify points that are duplicated?
  if (functionSpace.type() == "StructuredColumns") {
    if (grid.name().compare(0, 1, std::string{"L"}) == 0) {
      const atlas::functionspace::StructuredColumns fs(functionSpace);
      const atlas::StructuredGrid sgrid = fs.grid();
      const auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
      const auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
      for (atlas::idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
        for (atlas::idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
          atlas::idx_t jnode = fs.index(i, j);
          if (((view_j(jnode) == 1) || (view_j(jnode) == sgrid.ny())) && (view_i(jnode) != 1)) {
            ownedView(jnode, 0) = 0;
          }
        }
      }
    }
  }

  fieldset.add(owned);
}

// -----------------------------------------------------------------------------

}  // namespace util
