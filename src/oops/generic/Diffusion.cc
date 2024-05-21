/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <algorithm>
#include <utility>

#include "atlas/mesh/actions/BuildEdges.h"

#include "oops/base/GeometryData.h"
#include "oops/generic/Diffusion.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {

// --------------------------------------------------------------------------------------

Diffusion::Diffusion(const GeometryData & geometryData)
  : geom_(geometryData)
{
  oops::Log::trace() << "Diffusion::Diffusion start" << std::endl;
  util::Timer timer("oops::Diffusion", "Diffusion");

  // to avoid floating point exceptions, this is the smallest grid length we'll consider
  const double MIN_LENGTH = 1e-6;

  // Get the mesh connectivity information we need. This only works with
  // NodeColumns at the moment.
  // TODO(Travis) add support for StructuredColumns
  const atlas::FunctionSpace & fs = geom_.functionSpace();
  if (fs.type() != "NodeColumns") {
    util::abor1_cpp("oops::Diffusion : only NodeColumns atlas function space is supported for now",
                    __FILE__, __LINE__);
  }
  atlas::Mesh mesh(atlas::functionspace::NodeColumns(fs).mesh());
  atlas::mesh::actions::build_edges(mesh);
  ASSERT(mesh.nodes().has_field("xyz"));
  ASSERT(mesh.cells().has_field("centre"));
  ASSERT(mesh.edges().size() > 0);

  // process the geometry. (Noting that the atlas mesh nodes are the center of
  // our model grid cell) For each pair of connecting nodes we need to calculate
  // the length of the edge connecting two nodes (easy) and the length of the
  // model grid box side that crosses this edge (not as easy).
  // ----------------------------------------------------------------------------------------------

  // get the fields we'll need later
  const auto ghost = atlas::array::make_view<int, 1>(mesh.nodes().ghost());
  const auto xyz = atlas::array::make_view<double, 2>(mesh.nodes().field("xyz"));
  const auto centers = atlas::array::make_view<double, 2>(mesh.cells().field("centre"));

  // get the edge/node/cell connectivity
  const auto & edge2node = mesh.edges().node_connectivity();
  const auto & edge2cell = mesh.edges().cell_connectivity();

  // calculate the grid parameters we'll need later for diffusion.
  edgeGeom_.reserve(mesh.edges().size());
  for (atlas::idx_t i = 0; i < mesh.edges().size(); i++) {
    // get the node indexes
    ASSERT(edge2node.cols(i) == 2);
    const auto nodeA = edge2node(i, 0);
    const auto nodeB = edge2node(i, 1);
    ASSERT(nodeA != nodeB);

    // If both nodes are in the halo, don't bother adding this edge to the vector
    if (ghost(nodeA) && ghost(nodeB)) continue;

    // create edge in vector. Make sure the lowest value index is first for
    // reasons that I might care about later, maybe.
    EdgeGeom &edgeGeom = edgeGeom_.emplace_back();
    edgeGeom.nodeA = nodeA;
    edgeGeom.nodeB = nodeB;
    if (edgeGeom.nodeA > edgeGeom.nodeB) {
      std::swap(edgeGeom.nodeA, edgeGeom.nodeB);
    }

    // calculate the length of the mesh edge (i.e. length between two model grid cell centers)
    const auto & pointA = atlas::Point3(xyz(nodeA, 0), xyz(nodeA, 1), xyz(nodeA, 2));
    const auto & pointB = atlas::Point3(xyz(nodeB, 0), xyz(nodeB, 1), xyz(nodeB, 2));
    edgeGeom.edgeLength = atlas::Point3::distance(pointA, pointB);

    // get atlas mesh cell centers, and estimate length of original model grid cell edge
    // that passes through this atlas edge.
    ASSERT(edge2cell.cols(i) == 2);
    const atlas::idx_t cellA = edge2cell(i, 0);
    const atlas::idx_t cellB = edge2cell(i, 1);

    // NOTE atlas is returning a cell index of -1 if there is only 1 cell. Also,
    // at one point with some compilers I was getting very large cell indexes
    // for invalid cells, not sure if that is still a problem
    if (cellA < 0 || cellB < 0 ||
        cellA >= mesh.cells().size() || cellB >= mesh.cells().size()) {
      // There is no second cell center, so just make up a reasonable value
      edgeGeom.aspectRatio = 1.0;
    } else {
      const auto & centerA = atlas::Point3(centers(cellA, 0),
                                           centers(cellA, 1),
                                           centers(cellA, 2));
      const auto & centerB = atlas::Point3(centers(cellB, 0),
                                           centers(cellB, 1),
                                           centers(cellB, 2));
      const double center2centerLen = atlas::Point3::distance(centerA, centerB);
      edgeGeom.aspectRatio = edgeGeom.edgeLength < MIN_LENGTH ? 0.0 :
                             center2centerLen / edgeGeom.edgeLength;
    }
  }

  // save other grid based constants that are used later. e.g. inv_area (1/area)
  inv_area_ = fs.createField<double>();
  auto v_inv_area = atlas::array::make_view<double, 1>(inv_area_);
  if (geom_.fieldSet().has("area")) {
    // yay, the model interface gave us an area.
    auto v_area = atlas::array::make_view<double, 2>(geom_.fieldSet().field("area"));
    for (atlas::idx_t i = 0; i < fs.size(); i++) {
      v_inv_area(i) = v_area(i, 0) < MIN_LENGTH ? 0.0 : 1.0 / v_area(i, 0);
    }
  } else {
    // Ugh, the model interface did NOT give us any precomputed area field. Warn
    // user that diffusion accuracy might be compromised.
    oops::Log::warning() << "WARNING - oops::Diffusion needs the geometry's area but none is"
                            " given. Crudely estimating area. Please provide"
                            " an 'area' field in your geometry." << std::endl;

    // just get a rough estimate of the the area, it doesn't have to be exact.
    // For each edge, add its length to the two nodes, take the square of the
    // average of these distances for all the nodes
    auto count = fs.createField<int>();
    auto v_count = atlas::array::make_view<int, 1>(count);
    v_count.assign(0);
    v_inv_area.assign(0.0);
    for (const auto & edge : edgeGeom_) {
      v_count(edge.nodeA) += 1;
      v_count(edge.nodeB) += 1;
      v_inv_area(edge.nodeA) += edge.edgeLength;
      v_inv_area(edge.nodeB) += edge.edgeLength;
    }
    for (atlas::idx_t i = 0; i < fs.size(); i++) {
      if (v_count(i) == 0) continue;
      double v = v_inv_area(i) / v_count(i);
      v_inv_area(i) = 1.0 / (v * v);
    }
    inv_area_.set_dirty();
    inv_area_.haloExchange();

    // TODO(Travis) if we stop having models provide area, look into having a
    // more accurate area estimate.
  }

  oops::Log::trace() << "Diffusion::Diffusion end" << std::endl;
}

// --------------------------------------------------------------------------------------

void Diffusion::setParameters(const atlas::FieldSet & parameters) {
  oops::Log::trace() << "Diffusion::setParameters start" << std::endl;
  util::Timer timer("oops::Diffusion", "setParameters");

  const std::string HZ_SCALES = "hzScales";
  const std::string VT_SCALES = "vtScales";

  const atlas::FunctionSpace & fs = geom_.functionSpace();
  parameters.haloExchange();

  // sanity check
  if (!parameters.has(HZ_SCALES) && !parameters.has(VT_SCALES)) {
    util::abor1_cpp("Diffusion::setScales() neither 'hzScales' nor 'vtScales' was given.",
      __FILE__, __LINE__);
  }

  //-------------------------------------------------------------------------------------
  // Horizontal diffusion parameters (in units of m)
  //-------------------------------------------------------------------------------------
  if (parameters.has(HZ_SCALES)) {
    // calculate the actual min number of hz iterations, and the diffusion coefficients (khdt_)
    double minItr = 0;
    const auto v_hzScales = atlas::array::make_view<double, 2>(parameters[HZ_SCALES]);
    khdtLevels_ = v_hzScales.shape(1);
    khdt_.resize(edgeGeom_.size(), std::vector<double>(khdtLevels_));

    for (size_t e = 0; e < khdt_.size(); e++) {
      for (int level = 0; level < khdtLevels_; level++) {
        if (v_hzScales(edgeGeom_[e].nodeA, level) == 0.0 ||
            v_hzScales(edgeGeom_[e].nodeB, level) == 0.0) {
          // one of the nodes is masked out, dont diffuse with it
          khdt_[e][level] = 0.0;
        } else {
          // calculate diffusion coefficient (not taking into account the number of
          // iterations.. yet)
          const double s = (v_hzScales(edgeGeom_[e].nodeA, level) +
                            v_hzScales(edgeGeom_[e].nodeB, level)) / 2.0;
          khdt_[e][level] = s * s;

          // calculate the minimum number of iterations needed to be computationally stable
          // on this PE
          const double el = edgeGeom_[e].edgeLength;
          if (el <= 0) {
            util::abor1_cpp("oops::Diffusion will not work with grids with degenerate points"
                            " that are not masked out", __FILE__, __LINE__);
          }
          minItr = std::max(2.0 * (s*s) / (el*el), minItr);
        }
      }
    }
    niterHz_ = std::ceil(minItr/2) * 2;  // make sure number of iterations is even
    // get the global min number of iterations
    geom_.comm().allReduceInPlace(niterHz_, eckit::mpi::Operation::MAX);

    // TODO(Travis) do some error checking to make sure niterHz_ is not too big
    oops::Log::info() << "Diffusion: horizontal iterations: " << niterHz_ << std::endl;

    // adjust the above calculated diffusion coefficients by the final number of
    // iterations
    for (size_t e = 0; e < khdt_.size(); e++) {
      for (int level = 0; level < khdtLevels_; level++) {
        khdt_[e][level] *= 1.0 / (2.0 * niterHz_);
      }
    }
  }

  //-------------------------------------------------------------------------------------
  // Vertical diffusion parameters (in units of # of levels)
  //-------------------------------------------------------------------------------------
  if (parameters.has(VT_SCALES)) {
    // calculate the min number of vertical iterations, and the diffusion coefficients
    double minItr = 0;
    const auto v_vtScales = atlas::array::make_view<double, 2> (parameters[VT_SCALES]);
    kvdt_ = fs.createField<double>(atlas::option::levels(v_vtScales.shape(1)-1));
    auto v_kvdt = atlas::array::make_view<double, 2>(kvdt_);
    for (atlas::idx_t i = 0; i < kvdt_.shape(0); i++) {
      for (atlas::idx_t level = 0; level < kvdt_.shape(1); level++) {
        if (v_vtScales(i, level) == 0.0 || v_vtScales(i, level+1) == 0.0) {
          // one of the nodes is masked out, don't diffuse with it.
          v_kvdt(i, level) = 0.0;
        } else {
          // calculate the diffusion coefficient (not taking into account the number of
          // iterations.. yet)
          double s = (v_vtScales(i, level) + v_vtScales(i, level+1)) / 2;
          v_kvdt(i, level) = s * s;

          // calculate the minimum number of iterations needed to be computationally stable
          // on this PE
          minItr = std::max(2.0 * s * s, minItr);
        }
      }
    }
    niterVt_ = std::ceil(minItr/2) * 2;  // make sure number of iterations is even
    // get the global min number of iterations
    geom_.comm().allReduceInPlace(niterVt_, eckit::mpi::Operation::MAX);

    // TODO(Travis) do some error checking to make sure niterVt_ is not too big
    oops::Log::info() << "Diffusion: vertical iterations: " << niterVt_ << std::endl;

    // adjust the above calculated diffusion coefficients by the final number of
    // iterations
    for (atlas::idx_t i = 0; i < kvdt_.shape(0); i++) {
      for (atlas::idx_t lvl = 0; lvl < kvdt_.shape(1); lvl++) {
        v_kvdt(i, lvl) *= 1.0 / (2.0 * niterVt_);
      }
    }
  }
  oops::Log::trace() << "Diffusion::setScales end" << std::endl;
}

// --------------------------------------------------------------------------------------

void Diffusion::multiply(atlas::FieldSet &fset, Mode mode) const {
  oops::Log::trace() << "Diffusion::multiply start" << std::endl;
  util::Timer timer("oops::Diffusion", "multiply");

  // We do not have the true 3D diffusion operator implemented (yet?)
  // force user to use the split Hz/vt instead
  if (mode == Mode::Joint3D) {
    throw eckit::NotImplemented(
      "Diffusion:multiply does not have \"HZVT_3D\" mode implemented, use \"HZVT_2D_1D instead.",
      Here());
  }

  for (atlas::Field & field : fset) {
    bool doVt = (mode == Mode::VerticalOnly || mode == Mode::Split3D) && niterVt_ > 0;
    bool doHz = (mode == Mode::HorizontalOnly || mode == Mode::Split3D) && niterHz_ > 0;

    // apply half the vertical diffusion iterations
    if (doVt) multiplyVtTL(field);

    // then all the horizontal diffusion iterations
    if (doHz) {
     multiplyHzTL(field);
     multiplyHzTL(field);
    }

    // then the other half of the vertical diffusion iterations
    if (doVt) multiplyVtTL(field);
  }

  oops::Log::trace() << "Diffusion::multiply end" << std::endl;
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplySqrtTL(atlas::FieldSet &fset, Mode mode) const {
  oops::Log::trace() << "Diffusion::multiplySqrtTL start" << std::endl;
  util::Timer timer("oops::Diffusion", "multiplySqrtTL");

  // We do not have the true 3D diffusion operator implemented (yet?)
  // force user to use the split Hz/vt instead
  if (mode == Mode::Joint3D) {
    throw eckit::NotImplemented(
      "Diffusion:multiply does not have \"HZVT_3D\" mode implemented, use \"HZVT_2D_1D instead.",
      Here());
  }

  // for each field, diffuse!
  for (atlas::Field & field : fset) {
    bool doVt = (mode == Mode::VerticalOnly || mode == Mode::Split3D) && niterVt_ > 0;
    bool doHz = (mode == Mode::HorizontalOnly || mode == Mode::Split3D) && niterHz_ > 0;

    if (doHz) multiplyHzTL(field);
    if (doVt) multiplyVtTL(field);
  }
  oops::Log::trace() << "Diffusion::multiplySqrtTL end" << std::endl;
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplySqrtAD(atlas::FieldSet &fset, Mode mode) const {
  oops::Log::trace() << "Diffusion::multiplySqrtAD start" << std::endl;
  util::Timer timer("oops::Diffusion", "multiplySqrtAD");

  // We do not have the true 3D diffusion operator implemented (yet?)
  // force user to use the split Hz/vt instead
  if (mode == Mode::Joint3D) {
    throw eckit::NotImplemented(
      "Diffusion:multiply does not have \"HZVT_3D\" mode implemented, use \"HZVT_2D_1D instead.",
      Here());
  }

  // for each field, diffuse!
  for (atlas::Field & field : fset) {
    bool doVt = (mode == Mode::VerticalOnly || mode == Mode::Split3D) && niterVt_ > 0;
    bool doHz = (mode == Mode::HorizontalOnly || mode == Mode::Split3D) && niterHz_ > 0;

    if (doVt) multiplyVtAD(field);
    if (doHz) multiplyHzAD(field);
  }
  oops::Log::trace() << "Diffusion::multiplySqrtAD end" << std::endl;
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplyHzTL(atlas::Field & field) const {
  const atlas::FunctionSpace & fs = geom_.functionSpace();

  ASSERT(field.shape(0) == fs.size());
  ASSERT(khdtLevels_ == 1 || field.shape(1) <= khdtLevels_);

  const auto inv_area = atlas::array::make_view<double, 1>(inv_area_);
  auto fieldVal = atlas::array::make_view<double, 2>(field);
  std::vector<double> flux(edgeGeom_.size(), 0.0);

  for (atlas::idx_t itr = 0; itr < niterHz_/2; itr++) {
    field.haloExchange();

    for (atlas::idx_t level = 0; level < field.shape(1); level++) {
      // calculate diffusive flux at each edge. khdtLevels is the level to pull
      // khdt from, this code will work with either a 3D khdt field or a 2D khdt
      // field.
      const atlas::idx_t khdtLevel = std::min(level, khdtLevels_-1);
      for (size_t e = 0; e < edgeGeom_.size(); e++) {
        const double dv = fieldVal(edgeGeom_[e].nodeA, level) - fieldVal(edgeGeom_[e].nodeB, level);
        flux[e] = (edgeGeom_[e].aspectRatio * dv) * khdt_[e][khdtLevel];
      }

      // time-step the diffusion terms
      for (size_t e = 0; e < edgeGeom_.size(); e++) {
        fieldVal(edgeGeom_[e].nodeA, level) -= inv_area(edgeGeom_[e].nodeA) * flux[e];
        fieldVal(edgeGeom_[e].nodeB, level) += inv_area(edgeGeom_[e].nodeB) * flux[e];
      }
    }
    field.set_dirty(true);
  }
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplyHzAD(atlas::Field & field) const {
  const atlas::FunctionSpace & fs = geom_.functionSpace();

  ASSERT(field.shape(0) == fs.size());
  ASSERT(khdtLevels_ == 1 || field.shape(1) <= khdtLevels_);

  const auto & inv_area = atlas::array::make_view<double, 1>(inv_area_);
  const auto & ghost = atlas::array::make_view<int, 1>(fs.ghost());
  auto fieldVal = atlas::array::make_view<double, 2>(field);

  std::vector<double> flux(edgeGeom_.size(), 0.0);

  // init halo to zero
  for (atlas::idx_t i = 0; i < fs.size(); i++) {
    if (ghost(i)) {
      for (atlas::idx_t level = 0; level < field.shape(1); level++) {
        fieldVal(i, level) = 0;
      }
    }
  }

  for (int itr = 0; itr < niterHz_/2; itr++) {
    for (int level = 0; level < field.shape(1); level++) {
      // adjoint time-step the diffusion terms
      for (size_t e = 0; e < edgeGeom_.size(); e++) {
        flux[e] += inv_area(edgeGeom_[e].nodeB) * fieldVal(edgeGeom_[e].nodeB, level)
                  -inv_area(edgeGeom_[e].nodeA) * fieldVal(edgeGeom_[e].nodeA, level);
      }

      // Adjoint calculate diffusive flux at each edge. khdtLevels is the level
      // to pull khdt from, this code will work with either a 3D khdt field or a
      // 2D khdt field.
      const int khdtLevel = std::min(level, khdtLevels_-1);
      for (size_t e = 0; e < edgeGeom_.size(); e++) {
        fieldVal(edgeGeom_[e].nodeA, level) += edgeGeom_[e].aspectRatio
                                               * khdt_[e][khdtLevel] * flux[e];
        fieldVal(edgeGeom_[e].nodeB, level) -= edgeGeom_[e].aspectRatio
                                               * khdt_[e][khdtLevel] * flux[e];
        flux[e] = 0.0;
      }
    }

    field.adjointHaloExchange();
  }

  field.set_dirty(true);
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplyVtTL(atlas::Field & field) const {
  if (field.shape(1) <= 1) return;  // early exit for 2D fields

  // make sure input field is correct shape
  ASSERT(field.shape(0) == kvdt_.shape(0));
  ASSERT(field.shape(1) == kvdt_.shape(1)+1);

  field.haloExchange();
  auto v_field = atlas::array::make_view<double, 2>(field);
  const auto & v_kvdt = atlas::array::make_view<double, 2>(kvdt_);

  const int nz = field.shape(1);
  std::vector<double> flux(nz-1, 0.0);  // 1 less level than field
                                        // since this is flux between levels

  for (int itr = 0; itr < niterVt_ / 2; itr++) {
    for (atlas::idx_t i = 0; i < field.shape(0); i++) {
      // calculate diffusive flux
      for (size_t level = 0; level < flux.size(); level++) {
        flux[level] = v_kvdt(i, level) * (v_field(i, level) - v_field(i, level+1));
      }

      // time-step diffusion terms
      for (int level = 0; level < nz-1; level++) {
        v_field(i, level) -= flux[level];
        v_field(i, level+1) += flux[level];
      }
    }
  }
}

// --------------------------------------------------------------------------------------
// NOTE: the vertical adjoint code should produce identical answers compared
// with the TL code above. It's explicitly coded anyway just to be safe
void Diffusion::multiplyVtAD(atlas::Field & field) const {
  if (field.shape(1) <= 1) return;  // early exit for 2D fields

  // make sure input field is correct shape
  ASSERT(field.shape(0) == kvdt_.shape(0));
  ASSERT(field.shape(1) == kvdt_.shape(1)+1);

  field.haloExchange();
  auto v_field = atlas::array::make_view<double, 2>(field);
  const auto & v_kvdt = atlas::array::make_view<double, 2>(kvdt_);

  const int nz = field.shape(1);
  std::vector<double> flux(nz-1, 0.0);  // 1 less level than field
                                        // since this is flux between levels

  for (int itr = 0; itr < niterVt_ / 2; itr++) {
    for (atlas::idx_t i = 0; i < field.shape(0); i++) {
      // adjoint time-step diffusion terms
      for (int level = 0; level < nz-1; level++) {
        flux[level] += v_field(i, level+1) - v_field(i, level);
      }

      // adjoint diffusive flux
      for (size_t level = 0; level < flux.size(); level++) {
        v_field(i, level) += v_kvdt(i, level) * flux[level];
        v_field(i, level+1) -= v_kvdt(i, level) * flux[level];
        flux[level] = 0.0;
      }
    }
  }
}

// --------------------------------------------------------------------------------------

}  // namespace oops
