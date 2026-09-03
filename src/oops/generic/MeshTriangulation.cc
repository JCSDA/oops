/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/generic/MeshTriangulation.h"

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation/element/Triag3D.h"
#include "atlas/interpolation/method/Ray.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/util/Point.h"

#include "eckit/mpi/Comm.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Timer.h"

namespace detail {
using atlas::Point3;
bool pointOutsideCircumcircle(const Point3 & point, const Point3 & a,
                              const Point3 & b, const Point3 & c) {
  atlas::PointXYZ center = Point3::cross(a, b) + Point3::cross(b, c) + Point3::cross(c, a);
  const double normA = Point3::norm(a);
  const double normCenter = Point3::norm(center);

  // Protect against the degenerate case where points coincide and the triangle has zero area.
  // This can occur for grids with singularities. If the triangle has zero area, then the
  // circumcircle is not well defined anyway, so we arbitrarily choose to return false. By
  // returning early, we avoid a div-by-zero when re-scaling `center` below.
  // Note: 0.5*norm(center) gives the triangle area; take sqrt to go back to linear units.
  if (sqrt(0.5 * normCenter) < 1e-14 * normA) {
    return false;
  }

  // Scale center to lie on sphere
  center *= normA / normCenter;

  // Check if point is outside circumcircle by comparing chord distances from center
  const double chordPoint = Point3::distance(point, center);
  const double chordCircle = Point3::distance(a, center);
  return chordPoint > (1.0 + 1e-6) * chordCircle;  // allow small sloppiness
}

bool checkFirstTriangulationIsDelaunay(const Point3 & p0, const Point3 & p1,
                                       const Point3 & p2, const Point3 & p3) {
  // First triangulation has triangles (p0,p1,p2) and (p2,p3,p0)
  // The Delaunay condition is that no other nodes lie in a triangle's circumcircle, so we check
  // that p3 is NOT in circumcircle of (p0,p1,p2) AND that p1 is NOT in circumcircle of (p2,p3,p0).
  return pointOutsideCircumcircle(p3, p0, p1, p2) && pointOutsideCircumcircle(p1, p2, p3, p0);
}
}  // namespace detail

namespace oops {

// -----------------------------------------------------------------------------

MeshTriangulation::MeshTriangulation(const atlas::FunctionSpace & fspace,
                                     const eckit::mpi::Comm & comm):
  firstTriangulationOfQuadsIsDelaunay_(), earth_(atlas::util::Earth::radius()),
  localCellCenterTree_(earth_)
{
  // Set default communicator name
  eckit::mpi::setCommDefault(comm.name().c_str());

  // Error for configurations that don't support triangulation
  ASSERT(fspace);  // check was initialized
  ASSERT(fspace.type() != "Spectral");

  // Check for custom MPI partitions where some tasks handle zero points
  // This configuration could likely be supported after adding new logic to skip/handle work on
  // the zero-size meshes that arise, but we haven't done that yet as it's a bit of an edge case.
  int min_nb_cells = fspace.size();
  comm.allReduceInPlace(min_nb_cells, eckit::mpi::min());
  ASSERT_MSG(min_nb_cells > 0,
             "MeshTriangulation received FunctionSpace " + fspace.type()
             + " using a distribution where some MPI tasks own zero points.");

  // Initialize mesh and local cell-center tree
  setMeshAndTriangulation(fspace, comm);
  setLocalTree();

  // This block must come after the mesh_ has been fully initialized
  if (fspace.type() == "StructuredColumns") {
    is_atlas_structured_columns_ = true;
    const atlas::functionspace::StructuredColumns structuredcolumns(fspace);
    indexMapper_.initialize(mesh_, structuredcolumns);

    const atlas::RegularGrid rg(structuredcolumns.grid());
    if (rg) {
      // Save regular grid nx for pole correction
      regular_grid_nx_ = rg.nx();
    }

  } else if (fspace.type() == "NodeColumns") {
    const atlas::RegularGrid rg(mesh_.grid());
    if (rg) {
      // Save regular grid nx for pole correction
      regular_grid_nx_ = rg.nx();
    }
  }
}

// -----------------------------------------------------------------------------

bool MeshTriangulation::containingTriangleAndBarycentricCoords(
    const double lat, const double lon,
    std::array<int, 3> & indices, std::array<double, 3> & baryCoords) const {
  ASSERT(mesh_);
  ASSERT(mesh_.nodes().has_field("xyz"));
  ASSERT(!localCellCenterTree_.empty());
  util::Timer timer("oops::MeshTriangulation", "containingTriangleAndBarycentricCoords");

  const auto fillReturnArgsFromIntersect = [&](
      const atlas::interpolation::method::Intersect & intersect,
      const int indexA, const int indexB, const int indexC) {
    ASSERT(intersect);
    if (is_atlas_structured_columns_) {
      indices[0] = indexMapper_(indexA);
      indices[1] = indexMapper_(indexB);
      indices[2] = indexMapper_(indexC);
    } else {
      indices[0] = indexA;
      indices[1] = indexB;
      indices[2] = indexC;
    }
    baryCoords[0] = 1.0 - intersect.u - intersect.v;
    baryCoords[1] = intersect.u;
    baryCoords[2] = intersect.v;
    // The atlas coordinates u,v are in [0,1], but can still have roundoff-level negative 1-u-v
    if (baryCoords[0] < 0.0) {
      ASSERT(fabs(baryCoords[0]) < 1e-14);  // negative but larger than roundoff is a bug
      baryCoords[0] = 0.0;
    }
  };

  const auto & connectivity = mesh_.cells().node_connectivity();
  const auto & xyz = atlas::array::make_view<double, 2>(mesh_.nodes().field("xyz"));

  const auto makePoint3 = [&](const int localindex) -> atlas::Point3 {
    return atlas::Point3(xyz(localindex, 0), xyz(localindex, 1), xyz(localindex, 2));
  };

  const auto checkPointInSphericalTriangle = [&](const atlas::Point3 & point,
      const int cell, const int nodeA, const int nodeB, const int nodeC) -> bool {
    const int indexA = connectivity(cell, nodeA);
    const int indexB = connectivity(cell, nodeB);
    const int indexC = connectivity(cell, nodeC);

    // Fail early if any vertex of this triangle would fail to be handled by the index remapper.
    // This indicates that we're deep enough into the Mesh's halo that the FunctionSpace does not
    // have a matching grid point, so no meaningful stencil can be generated. This scenario can
    // arise when using atlas StructuredColumns distributed over a small number of MPI ranks: in
    // this configuration halos can overlap owned points across the lon=0 meridian, leading to
    // multiple triangles containing the target point, but some will be deep in the Mesh halo so
    // must be discarded.
    if (is_atlas_structured_columns_) {
      const auto missing = util::missingValue<atlas::idx_t>();
      if (indexMapper_(indexA) == missing
          || indexMapper_(indexB) == missing
          || indexMapper_(indexC) == missing) {
        return false;
      }
    }

    const atlas::Point3 a = makePoint3(indexA);
    const atlas::Point3 b = makePoint3(indexB);
    const atlas::Point3 c = makePoint3(indexC);
    const atlas::interpolation::element::Triag3D tri(a, b, c);
    const double sqrtArea = sqrt(tri.area());

    // Protect against the degenerate case where points coincide and the triangle has zero area.
    // This can occur for grids with singularities. If the triangle has zero area, we return early
    // and make sure it does not contain the target point.
    // Note: atlas computes XYZ mesh points using Earth geometry, so the 3D triangle lives on the
    // Earth's surface. We take the sqrt to go back to linear units, and compare to Earth's radius.
    if (sqrtArea < 1e-14 * earth_.radius()) {
      return false;
    }

    const atlas::interpolation::method::Ray ray(point);
    const double edgeEpsilon = 1e-15 * sqrtArea;
    const auto intersect = tri.intersects(ray, edgeEpsilon);

    if (intersect) {
      fillReturnArgsFromIntersect(intersect, indexA, indexB, indexC);
      return true;
    } else {
      return false;
    }
  };

  // The number of cells to check depends on how the search pattern interacts with the mesh
  const int nb_cells_to_check = [&]() {
    // Default case: check the 8 closest cells, following the atlas unstructured interpolator
    int nb_to_check = 8;
    // Special case: in a regular grid, we expect the convergence of grid lines at the pole will
    // lead to highly-deformed cells (quads or tris). This makes it necessary to search through more
    // cells than the default case to guarantee finding the cell containing the target point, which
    // in turn makes our proximity-based search inefficient. If this expanded search does become a
    // bottleneck, an alternative search pattern in polar regions should be considered; perhaps
    // using cell-to-cell connectivity to search through neighbors of neighbors.
    if (regular_grid_nx_ > 0) {
      // Scale number of cells by the aspect ratio of the grid cells, which is approximately given
      // by the compression of the const-longitude lines towards the pole:
      const double max_to_check = 2.0 * regular_grid_nx_;  // arbitrary max: 2 bands of cells
      const double deg2rad = M_PI / 180.0;
      const double compression = 1.0 / std::max(std::abs(std::cos(deg2rad * lat)), 1e-14);
      // Apply scaling factor and check against max; we do this as floating-point math, because if
      // the target lat is close to a pole, then the compression will be huge and integers overflow.
      const double nb_scaled_and_bounded = std::min(nb_to_check * compression, max_to_check);
      nb_to_check = static_cast<int>(nb_scaled_and_bounded);
    }
    return std::min(nb_to_check, static_cast<int>(localCellCenterTree_.size()));
  }();

  // Sort the list of cells returned by KD-tree. The list is already sorted by distance,
  // sort only equidistant points by payload (cell index).
  const auto sortOnlyTies = [](auto & list) {
    auto first = list.begin();
    while (first != list.end()) {
      auto rangeEnd = std::find_if(first, list.end(),
                [d = first->distance()](const auto & x) { return x.distance() != d; });
      if (std::distance(first, rangeEnd) > 1) {
        std::sort(first, rangeEnd, [](const auto & a, const auto & b) {
                                   return a.payload() < b.payload();
                                   });
      }
      first = rangeEnd;
    }
  };

  // Find cell that contains target point
  atlas::PointLonLat pll(lon, lat);
  pll.normalise();
  atlas::Point3 p;
  earth_.lonlat2xyz(pll, p);

  bool success = false;

  auto list = localCellCenterTree_.closestPoints(p, nb_cells_to_check);
  // The list is already sorted by distance, now sort only equidistant points by
  // payload (cell index). This ensures that the order is deterministic and results are
  // reproducible with different MPI layouts.
  sortOnlyTies(list);
  for (const auto & item : list) {
    const int cell = item.payload();
    const int nb_cols = connectivity.cols(cell);

    if (nb_cols == 3) {
      success = checkPointInSphericalTriangle(p, cell, 0, 1, 2);
      if (success) { break; }
    } else {
      ASSERT(!firstTriangulationOfQuadsIsDelaunay_.empty());
      if (firstTriangulationOfQuadsIsDelaunay_[cell]) {
        // triangle (p0,p1,p2)
        success = checkPointInSphericalTriangle(p, cell, 0, 1, 2);
        if (success) { break; }
        // triangle (p2,p3,p0)
        success = checkPointInSphericalTriangle(p, cell, 2, 3, 0);
        if (success) { break; }
      } else {
        // triangle (p3,p0,p1)
        success = checkPointInSphericalTriangle(p, cell, 3, 0, 1);
        if (success) { break; }
        // triangle (p1,p2,p3)
        success = checkPointInSphericalTriangle(p, cell, 1, 2, 3);
        if (success) { break; }
      }
    }
  }

  // The most likely explanation for failing to locate the target point is that it lies outside of
  // a regional grid. It is also possible (but unlikely) that the triangle containing the target
  // point is not in the group of triangles checked, i.e., is not one of the 8 closest triangles.
  // For now, we return false and let the client handle the failure to locate the target point.
  return success;
}

// -----------------------------------------------------------------------------

void MeshTriangulation::setMeshAndTriangulation(const atlas::FunctionSpace & fspace,
                                                const eckit::mpi::Comm & comm) {
  ASSERT(!mesh_);
  ASSERT(firstTriangulationOfQuadsIsDelaunay_.empty());
  util::Timer timer("oops::MeshTriangulationache", "setMeshAndTriangulation");

  // Setup mesh
  if (fspace.type() == "NodeColumns") {
    const atlas::functionspace::NodeColumns nodecolumns(fspace);
    mesh_ = nodecolumns.mesh();
  } else if (fspace.type() == "StructuredColumns") {
    const atlas::functionspace::StructuredColumns structuredcolumns(fspace);
    const atlas::StructuredGrid & grid = structuredcolumns.grid();
    if (fspace.distribution() == "custom") {
      // Gather global partition field on root processor
      atlas::Field globalPartition = fspace.createField<int>(
        atlas::option::name("partition") | atlas::option::global());
      fspace.gather(fspace.partition(), globalPartition);

      // Transform to a global partition vector
      std::vector<int> partition(grid.size());
      if (comm.rank() == 0) {
        ASSERT(grid.size() == static_cast<int>(globalPartition.size()));
        const auto globalPartitionView = atlas::array::make_view<int, 1>(globalPartition);
        for (atlas::idx_t jj = 0; jj < grid.size(); ++jj) {
          partition[jj] = globalPartitionView(jj);
        }
      }

      // Broadcast global partition vector
      comm.broadcast(partition, 0);

      // Create distribution and mesh
      atlas::grid::Distribution distribution;
      util::setupStructuredMeshWithCustomPartition(comm, grid, partition, distribution, mesh_);
    } else {
      // Create mesh
      const atlas::StructuredMeshGenerator gen(grid.meshgenerator());
      mesh_ = gen(grid, atlas::grid::Partitioner(fspace.distribution()));
    }

    // At this point, we have a mesh from the StructuredMeshGenerator that doesn't include a halo.
    // We add a halo via actions::build_halo, but BEWARE one critical caveat: the mesh halo is
    // structured DIFFERENTLY than the halo in the StructuredColumns FunctionSpace. In other words:
    // `fspace.lonlat()` will be a different set of points from `mesh.nodes().lonlat()` -- the same
    // owned points in the same order, but (in general) a different set of ghost points in a
    // different order. Therefore, one CANNOT use a mesh-based computation to determine an index
    // into the FunctionSpace/FieldSet, without first creating a mapping between the two halos.
    // See https://github.com/JCSDA-internal/oops/issues/2621
    if (grid.domain().global()) {
      // Atlas global structured grids, by default, are not periodic in longitude. We call the
      // action build_periodic_boundaries BEFORE build_halo to ensure the halos will cross the
      // lon=0 meridian.
      atlas::mesh::actions::build_nodes_parallel_fields(mesh_);
      atlas::mesh::actions::build_periodic_boundaries(mesh_);
    }
    atlas::mesh::actions::build_halo(mesh_, 1);

    // Reset default communicator name
    eckit::mpi::setCommDefault(comm.name().c_str());
  } else {
    ABORT(fspace.type() + " function space not supported yet");
  }

  // Add XYZ field to mesh, because these 3d coordinates are used in interpolation setup
  atlas::mesh::actions::BuildXYZField()(mesh_);

  const size_t nb_cells = mesh_.cells().size();
  const auto & connectivity = mesh_.cells().node_connectivity();
  const auto & xyz = atlas::array::make_view<double, 2>(mesh_.nodes().field("xyz"));

  const auto makePoint3 = [&](const int cell, const int node) -> atlas::Point3 {
    const int localindex = connectivity(cell, node);
    return atlas::Point3(xyz(localindex, 0), xyz(localindex, 1), xyz(localindex, 2));
  };

  // If model has quad cells, triangulate them.
  //
  // We do this split to satisfy the Delaunay condition locally on the quad, giving the better of
  // the two possible splitting diagonals. Note this may not produce a globally-optimal Delaunay
  // triangulation, depending on how the model quads are set up.
  //
  // first triangulation  => split into triangles (p0,p1,p2) and (p2,p3,p0)
  // second triangulation => split into triangles (p3,p0,p1) and (p1,p2,p3)
  for (size_t i = 0; i < nb_cells; ++i) {
    const int nb_cols = connectivity.cols(i);
    ASSERT(nb_cols == 3 || nb_cols == 4);

    if (nb_cols == 4) {
      // Allocate triangulation data on first need
      if (firstTriangulationOfQuadsIsDelaunay_.empty()) {
        firstTriangulationOfQuadsIsDelaunay_.assign(nb_cells, false);
      }

      const atlas::Point3 p0 = makePoint3(i, 0);
      const atlas::Point3 p1 = makePoint3(i, 1);
      const atlas::Point3 p2 = makePoint3(i, 2);
      const atlas::Point3 p3 = makePoint3(i, 3);
      firstTriangulationOfQuadsIsDelaunay_[i] =
          detail::checkFirstTriangulationIsDelaunay(p0, p1, p2, p3);
    }
  }
}

// -----------------------------------------------------------------------------

void MeshTriangulation::setLocalTree() {
  ASSERT(mesh_);
  ASSERT(mesh_.nodes().has_field("xyz"));
  ASSERT(localCellCenterTree_.empty());
  util::Timer timer("oops::MeshTriangulation", "setLocalTree");

  const size_t nb_cells = mesh_.cells().size();

  const atlas::Field & centersField = atlas::mesh::actions::BuildCellCentres()(mesh_);
  const auto centersView = atlas::array::make_view<double, 2>(centersField);
  std::vector<atlas::Point3> centers(nb_cells);
  std::vector<int> indices(nb_cells);
  for (size_t i = 0; i < nb_cells; ++i) {
    centers[i] = atlas::Point3(centersView(i, 0), centersView(i, 1), centersView(i, 2));
    indices[i] = i;
  }
  localCellCenterTree_.build(centers, indices);
}

// -----------------------------------------------------------------------------

namespace {

/// Flyweight factory to manage creation and retrieval of MeshTriangulation objects
class MeshTriangulationFactory {
 public:
  static const MeshTriangulation & get(const atlas::FunctionSpace & fspace,
                                                   const eckit::mpi::Comm & comm) {
    const std::string sep = "_";
    const std::string key = comm.name() + sep + util::getGridUid(fspace) + sep
                            + util::getLonLatHash(fspace);

    auto & map = instances();
    const auto it = map.find(key);
    if (it != map.end()) {
      return *(it->second);
    }
    const auto inserted = map.emplace(key,
        std::make_unique<MeshTriangulation>(fspace, comm));
    Log::debug() << "MeshTriangulationFactory: created new entry for key '"
                 << key << "'." << std::endl;
    return *(inserted.first->second);
  }

 private:
  MeshTriangulationFactory() {}
  MeshTriangulationFactory(const MeshTriangulationFactory &) = delete;
  MeshTriangulationFactory & operator=(const MeshTriangulationFactory &) = delete;

  // The map is intentionally heap-allocated and never deleted, leaking its contents to work
  // around a static-destruction-order crash: an atlas::Mesh destroyed after oops::Run has
  // finalised MPI and the atlas runtime corrupts the heap. A cleaner approach would register a
  // clear() with a library-level finalisation callback registry invoked before MPI finalisation.
  static std::unordered_map<std::string, std::unique_ptr<MeshTriangulation>> &
  instances() {
    static auto * instances =
        new std::unordered_map<std::string, std::unique_ptr<MeshTriangulation>>();
    return *instances;
  }
};

}  // namespace

const MeshTriangulation & getMeshTriangulation(
    const atlas::FunctionSpace & fspace, const eckit::mpi::Comm & comm) {
  return MeshTriangulationFactory::get(fspace, comm);
}

// -----------------------------------------------------------------------------

}  // namespace oops
