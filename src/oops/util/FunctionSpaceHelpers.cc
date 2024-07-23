/*
 * (C) Copyright 2024 UCAR
 * (C) Crown Copyright 2024 Met Office
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
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/meshgenerator.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/missingValues.h"

namespace util {


atlas::idx_t getSizeOwned(const atlas::FunctionSpace & fspace) {
  atlas::idx_t size_owned;
  executeFunc(fspace, [&](const auto& fspace){size_owned = fspace.sizeOwned();});
  return size_owned;
}

// -----------------------------------------------------------------------------

void StructuredMeshToStructuredColumnsIndexMap::initialize(
    const atlas::Mesh & mesh, const atlas::functionspace::StructuredColumns & structuredcolumns) {
  oops::Log::trace() << "StructuredMeshToStructuredColumnsIndexMap::initialize start" << std::endl;

  const auto mesh_lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
  const auto mesh_ghost = atlas::array::make_view<int, 1>(mesh.nodes().ghost());

  // Sanity checks
  nb_mesh_owned_ = 0;
  nb_mesh_ghost_ = 0;
  for (int jmesh = 0; jmesh < mesh.nodes().size(); ++jmesh) {
    if (mesh_ghost(jmesh) == 0) {
      ++nb_mesh_owned_;
      // sanity check: no ghost points come before any owned points
      ASSERT(nb_mesh_ghost_ == 0);
    } else {
      ++nb_mesh_ghost_;
    }
  }
  ASSERT(nb_mesh_owned_ + nb_mesh_ghost_ == mesh.nodes().size());

  const auto fs_lonlat = atlas::array::make_view<double, 2>(structuredcolumns.lonlat());
  const auto fs_ghost = atlas::array::make_view<int, 1>(structuredcolumns.ghost());
  atlas::idx_t nb_fs_owned = 0;
  atlas::idx_t nb_fs_ghost = 0;
  for (int jfs = 0; jfs < fs_ghost.shape(0); ++jfs) {
    if (fs_ghost(jfs) == 0) {
      ++nb_fs_owned;
      // sanity check: no ghost points come before any owned points
      ASSERT(nb_fs_ghost == 0);
    } else {
      ++nb_fs_ghost;
    }
  }
  ASSERT(nb_fs_owned + nb_fs_ghost == fs_ghost.shape(0));
  // sanity check: Mesh and FunctionSpace have same number of owned points
  ASSERT(nb_mesh_owned_ == nb_fs_owned);

  // Resize the map, filling with the missingValue denoting that no mapping was established
  map_.resize(nb_mesh_ghost_, missingValue<atlas::idx_t>());

  // Handle two trivial edge cases:
  // 1. Mesh has no ghost points, so there's no mapping to do or data to initialize
  if (nb_mesh_ghost_ == 0) {
    valid_ = true;
    return;
  }
  // 2. FunctionSpace has no ghost points, so all mesh ghosts map onto missing; since that's the
  //    default initialization, there's no further work to do
  if (nb_fs_ghost == 0) {
    valid_ = true;
    return;
  }

  // Build search tree for halo points' coordinates
  const atlas::Geometry earth(atlas::util::Earth::radius());
  atlas::util::IndexKDTree2D ghostTree(earth);
  std::vector<double> tree_lons(nb_fs_ghost);
  std::vector<double> tree_lats(nb_fs_ghost);
  std::vector<atlas::idx_t> tree_indices(nb_fs_ghost);
  for (int jghost = 0; jghost < nb_fs_ghost; ++jghost) {
    const atlas::idx_t fs_index = nb_fs_owned + jghost;
    ASSERT(fs_ghost(fs_index) == 1);
    tree_lons[jghost] = fs_lonlat(fs_index, 0);
    tree_lats[jghost] = fs_lonlat(fs_index, 1);
    tree_indices[jghost] = fs_index;
  }
  ghostTree.build(tree_lons, tree_lats, tree_indices);

  // Compute indices by looking up in the tree
  // Optimization note: some indices can be remapped using atlas's methods grid.index2ij then
  // structuredcolumns.index. The trouble is that won't work for indices in the "external" halos
  // at the "edge" of the domain, as views in the canonical [0,360]x[-90,90] coordinate patch.
  // We use the tree-based search because it's simple and handles all cases. A minor optimization
  // might be to use exact mapping for internal halos and tree-based search for "external" halos...
  for (atlas::idx_t jghost = 0; jghost < nb_mesh_ghost_; ++jghost) {
    const atlas::idx_t mesh_index = nb_mesh_owned_ + jghost;
    ASSERT(mesh_ghost(mesh_index) == 1);
    const atlas::PointLonLat pll(mesh_lonlat(mesh_index, 0), mesh_lonlat(mesh_index, 1));
    const auto item = ghostTree.closestPoint(pll);
    if (item.distance() < 1.e-9) {
      map_[jghost] = item.payload();
    }
  }

  valid_ = true;
  oops::Log::trace() << "StructuredMeshToStructuredColumnsIndexMap::initialize end" << std::endl;
}

// -----------------------------------------------------------------------------

atlas::idx_t StructuredMeshToStructuredColumnsIndexMap::operator()(
    const atlas::idx_t mesh_index) const {
  ASSERT(valid_);
  ASSERT(mesh_index >= 0 && mesh_index < nb_mesh_owned_ + nb_mesh_ghost_);
  if (mesh_index < nb_mesh_owned_) {
    return mesh_index;
  } else {
    const atlas::idx_t candidate = map_[mesh_index - nb_mesh_owned_];
    return candidate;
  }
}

// -----------------------------------------------------------------------------

void setupFunctionSpace(const eckit::mpi::Comm & comm,
    const eckit::Configuration & config,
    atlas::Grid & grid,
    atlas::grid::Partitioner & partitioner,
    atlas::Mesh & mesh,
    atlas::FunctionSpace & functionSpace,
    atlas::FieldSet & fieldset) {
  oops::Log::trace() << "setupFunctionSpace starting" << std::endl;
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

  const int halo = config.getInt("halo", 0);

  // Set up FunctionSpace (and Mesh)
  const std::string functionSpaceName = config.getString("function space");
  if (functionSpaceName == "StructuredColumns") {
    ASSERT(gridType != "unstructured");

    if (noPointOnLastTask && (comm.size() > 1)) {
      // Empirically, the atlas calls used here to create a custom distribution do NOT work when
      // requesting halos -- the function space constructor segfaults. If it becomes necessary to
      // use halos with the custom distribution, more investigation of the atlas interface may be
      // needed.
      if (halo > 0) {
        throw eckit::BadParameter(
            "Setting up a FunctionSpace with option `no point on last task` and with halos"
            " is currently not supported",
            Here());
      }

      // Create distribution from partitioner
      std::vector<int> partition(grid.size());
      partitioner.partition(grid, partition.data());

      // Create distribution and mesh
      atlas::grid::Distribution distribution;
      setupStructuredMeshWithCustomPartition(comm, grid, partition, distribution, mesh);

      // Create functionspace from distribution
      functionSpace = atlas::functionspace::StructuredColumns(grid, distribution,
                                                              atlas::option::halo(halo));
    } else {
      // Create mesh from partitioner
      mesh = atlas::MeshGenerator("structured").generate(grid, partitioner);

      // Create functionspace from partitioner
      functionSpace = atlas::functionspace::StructuredColumns(grid, partitioner,
                                                              atlas::option::halo(halo));

      // Using atlas::mpi::Scope in the call to atlas::functionspace::StructuredColumns
      // may have reverted the default communicator to the world communicator.
      // We set back the default communicator to `comm` to fix this.
      eckit::mpi::setCommDefault(comm.name().c_str());
    }

    // At this point, we have a mesh from the StructuredMeshGenerator that doesn't include a halo.
    // We add a halo via actions::build_halo, but BEWARE one critical caveat: the mesh halo is
    // structured DIFFERENTLY than the halo in the StructuredColumns FunctionSpace. In other words:
    // `fspace_.lonlat()` will be a different set of points from `mesh.nodes().lonlat()` -- the same
    // owned points in the same order, but (in general) a different set of ghost points in a
    // different order. Therefore, one CANNOT use a mesh-based computation to determine an index
    // into the FunctionSpace/FieldSet, without first creating a mapping between the two halos.
    // See https://github.com/JCSDA-internal/oops/issues/2621
    // To reduce the risk of incorrectly using the halos in the case of a structured mesh, we keep
    // the mesh halo-free for now.

    // Bugfix for regional grids
    // It seems that the content of lonlat for a regional function space is actually the xy
    // coordinates. The routine to compute distances on the sphere was complaining about
    // impossible lon/lat values...
    if (gridType == "regional") {
      auto lonlat = atlas::array::make_view<double, 2>(functionSpace.lonlat());
      double lonlatPoint[] = {0, 0};
      const atlas::functionspace::StructuredColumns fs(functionSpace);
      const atlas::StructuredGrid & grid = fs.grid();
      const auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
      const auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
      for (int jj = 0; jj < fs.size(); ++jj) {
        grid.lonlat(view_i(jj)-1, view_j(jj)-1, lonlatPoint);
        lonlat(jj, 0) = lonlatPoint[0];
        lonlat(jj, 1) = lonlatPoint[1];
      }
    }
  } else if (functionSpaceName == "NodeColumns") {
    if (grid.name().compare(0, 2, std::string{"CS"}) == 0) {
      // NodeColumns from a CubedSphere grid/mesh
      mesh = atlas::MeshGenerator("cubedsphere_dual",
                                  atlas::option::halo(halo)).generate(grid, partitioner);
      functionSpace = atlas::functionspace::CubedSphereNodeColumns(mesh);
    } else {
      // NodeColumns from an unstructured grid after triangulation
      // Regular or Structured grids could be supported with extra code
      mesh = atlas::MeshGenerator("delaunay").generate(grid, partitioner);
      atlas::mesh::actions::build_halo(mesh, halo);
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
  oops::Log::trace() << "setupFunctionSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

void setupStructuredMeshWithCustomPartition(const eckit::mpi::Comm & comm,
                                            const atlas::Grid & grid,
                                            const std::vector<int> & partition,
                                            atlas::grid::Distribution & distribution,
                                            atlas::Mesh & mesh) {
  // Create custom distribution
  std::vector<int> partitionCopy = partition;
  distribution = atlas::grid::Distribution(comm.size(), grid.size(), partitionCopy.data());

  // Count number of cells for each MPI task
  std::vector<size_t> nb_cells(comm.size(), 0);
  for (atlas::idx_t jj = 0; jj < grid.size(); ++jj) {
    ++nb_cells[partition[jj]];
  }

  // Get number of effective partitions (= number of tasks that own cells) and mapping from each
  // MPI rank to the effective partition number
  const int invalid_partition = comm.size() + 1;
  size_t nb_effective_partitions = 0;
  std::vector<size_t> effective_partition(comm.size(), invalid_partition);
  for (size_t jt = 0; jt < comm.size(); ++jt) {
    if (nb_cells[jt] > 0) {
      effective_partition[jt] = nb_effective_partitions;
      ++nb_effective_partitions;
    }
  }

  // Create mesh from distribution
  atlas::util::Config meshConfig(grid.meshgenerator());
  meshConfig.set("part", effective_partition[comm.rank()]);
  meshConfig.set("nb_parts", nb_effective_partitions);
  meshConfig.set("mpi_comm", comm.name());
  const atlas::StructuredMeshGenerator gen(meshConfig);
  mesh = gen(grid, distribution);
}

// -----------------------------------------------------------------------------

}  // namespace util
