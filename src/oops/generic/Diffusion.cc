/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <algorithm>
#include <set>
#include <utility>

#include "atlas/mesh/actions/BuildEdges.h"

#include "oops/base/GeometryData.h"
#include "oops/generic/Diffusion.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {

// --------------------------------------------------------------------------------------

// derived grid geometry. This struct is initialized at Diffusion constructor time,
// is independent of the desired diffusion scales, and can be copied in the
// Diffusion copy constructor. The way it is initialized depends on the type of
// atlas function space being used.
struct Diffusion::DerivedGeom {
  atlas::Field inv_area;  // Inverse of grid area
  struct EdgeGeom {
    size_t nodeA, nodeB;  // The two atlas nodes that this edge connects.
    double edgeLength;  // length between 2 atlas mesh nodes (i.e. between two model grid centers)
    double aspectRatio;  // edgeLength divided by  length of the perpendicular atlas cell centers
  };
  std::vector<EdgeGeom> edgeGeom;
};

// --------------------------------------------------------------------------------------

// Create the derived diffusion geometry from an atlas NodeColumn function
// space. The Mesh connectivity from NodeColumn is used to determine the derived
// geometry parameters
std::unique_ptr<Diffusion::DerivedGeom> calculateDerivedGeom_NodeColumns(
    const atlas::functionspace::NodeColumns & fs, const oops::GeometryData & geom) {
  std::unique_ptr<Diffusion::DerivedGeom> derivedGeom = std::make_unique<Diffusion::DerivedGeom>();

  // to avoid floating point exceptions, this is the smallest grid length we'll consider
  const double MIN_LENGTH = 1e-6;

  // make sure the mesh has what we need
  atlas::Mesh mesh(fs.mesh());
  try {
    atlas::mesh::actions::build_edges(mesh);
  } catch (...) {
    util::abor1_cpp("ERROR: oops::Diffusion() unable to build mesh edges. Make sure you"
                    " are running with > 1 PE, until oops is updated to use atlas version"
                    " >= 0.37", __FILE__, __LINE__);
  }
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
  derivedGeom->edgeGeom.reserve(mesh.edges().size());
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
    Diffusion::DerivedGeom::EdgeGeom &edgeGeom = derivedGeom->edgeGeom.emplace_back();
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
  derivedGeom->inv_area = fs.createField<double>();
  auto v_inv_area = atlas::array::make_view<double, 1>(derivedGeom->inv_area);
  if (geom.fieldSet().has("area")) {
    // yay, the model interface gave us an area.
    auto v_area = atlas::array::make_view<double, 2>(geom.fieldSet().field("area"));
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
    for (const auto & edge : derivedGeom->edgeGeom) {
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
    derivedGeom->inv_area.set_dirty();
    derivedGeom->inv_area.haloExchange();

    // TODO(Travis) if we stop having models provide area, look into having a
    // more accurate area estimate.
  }

  return derivedGeom;
}

// --------------------------------------------------------------------------------------

// Create the derived diffusion geometry from an atlas StructuredColumns function
// space. Note that we can't use Mesh in this case.
std::unique_ptr<Diffusion::DerivedGeom> calculateDerivedGeom_StructuredColumns(
    const atlas::functionspace::StructuredColumns & fs) {
  std::unique_ptr<Diffusion::DerivedGeom> derivedGeom = std::make_unique<Diffusion::DerivedGeom>();

  // Process the geometry. To keep the remainder of the Diffusion code generic
  // with respect to the NodeColumns function space, the StructuredColumns is
  // converted into a similar list of edges connecting two grid boxes. For each
  // edge, we calculate the length of the edge, and an estimate of the  length of
  // the model grid box side that crosses this edge.
  // ----------------------------------------------------------------------------------------------

  // make sure we have at least 1 halo available, otherwise this all breaks
  ASSERT(fs.halo() >= 1);

  // fields we'll need for later
  const auto lonlat = atlas::array::make_view<double, 2>(fs.lonlat());

  // calculate inverse area and temporary dx and dy fields
  auto dx = fs.createField<double>();
  auto vDx = atlas::array::make_view<double, 1>(dx);
  auto dy = fs.createField<double>();
  auto vDy = atlas::array::make_view<double, 1>(dy);
  derivedGeom->inv_area = fs.createField<double>();
  auto v_inv_area = atlas::array::make_view<double, 1>(derivedGeom->inv_area);
  eckit::geometry::SphereT<atlas::util::DatumIFS> earth;
  for (atlas::idx_t jj = fs.j_begin(); jj < fs.j_end(); jj++) {
    for (atlas::idx_t ii = fs.i_begin(jj); ii < fs.i_end(jj); ii++) {
      atlas::idx_t idx = fs.index(ii, jj);
      atlas::idx_t idx_im1 = fs.index(ii-1, jj);
      atlas::idx_t idx_ip1 = fs.index(ii+1, jj);
      atlas::idx_t idx_jm1 = fs.index(ii,   jj-1);
      atlas::idx_t idx_jp1 = fs.index(ii,   jj+1);

      // calculate dx / dy
      vDx(idx) = 0.5 * earth.distance(
        atlas::Point2(lonlat(idx_im1, 0), lonlat(idx_im1, 1)),
        atlas::Point2(lonlat(idx_ip1, 0), lonlat(idx_ip1, 1)));
      vDy(idx) = 0.5 * earth.distance(
        atlas::Point2(lonlat(idx_jm1, 0), lonlat(idx_jm1, 1)),
        atlas::Point2(lonlat(idx_jp1, 0), lonlat(idx_jp1, 1)));

      // calculate inverse area
      v_inv_area(idx) = 1.0 / (vDx(idx) * vDy(idx));
    }
  }
  dx.haloExchange();
  dy.haloExchange();
  derivedGeom->inv_area.haloExchange();

  // for each grid box (excluding halo), create information about the
  // connections with it's neighbors. Since we are iterating over grid points
  // (and not iterating over edges like in the NodeColumn version), we need to
  // keep track of which edges have already been counted, using "usedConnections"
  derivedGeom->edgeGeom.reserve(fs.sizeOwned()*3);  // we don't need that many, but just round up
  std::set<std::pair<atlas::idx_t, atlas::idx_t>> usedConnections;
  for (atlas::idx_t jj = fs.j_begin(); jj < fs.j_end(); jj++) {
    for (atlas::idx_t ii = fs.i_begin(jj); ii < fs.i_end(jj); ii++) {
      // get index of self and neighbors
      atlas::idx_t idx = fs.index(ii, jj);
      atlas::idx_t idx_im1 = fs.index(ii-1, jj);
      atlas::idx_t idx_ip1 = fs.index(ii+1, jj);
      atlas::idx_t idx_jm1 = fs.index(ii,   jj-1);
      atlas::idx_t idx_jp1 = fs.index(ii,   jj+1);

      // helper function to calculated the grid information we need between two
      // neighboring points
      auto addConnection = [&](auto idx1, auto idx2, auto &normField) {
        // find out if this connection has already been counted, if so, skip
        atlas::idx_t nodeA = idx1;
        atlas::idx_t nodeB = idx2;
        if (nodeA > nodeB) std::swap(nodeA, nodeB);
        std::pair<atlas::idx_t, atlas::idx_t> connection(nodeA, nodeB);
        if (usedConnections.find(connection) != usedConnections.end()) return;
        usedConnections.insert(connection);

        // create edge in vector
        Diffusion::DerivedGeom::EdgeGeom &edgeGeom = derivedGeom->edgeGeom.emplace_back();
        edgeGeom.nodeA = nodeA;
        edgeGeom.nodeB = nodeB;

        // calculate the distance between the two points
        atlas::Point2 pointA(lonlat(edgeGeom.nodeA, 0), lonlat(edgeGeom.nodeA, 1));
        atlas::Point2 pointB(lonlat(edgeGeom.nodeB, 0), lonlat(edgeGeom.nodeB, 1));
        edgeGeom.edgeLength = earth.distance(pointA, pointB);

        // calculate the aspect ratio. We already have the distance between the
        // points, estimate the length of the edge that crosses the edge between
        // those points
        double dNorm = 0.5 * (normField(idx1) + normField(idx2));
        edgeGeom.aspectRatio = dNorm / edgeGeom.edgeLength;
      };
      // Calculate grid information with all of our neighbors.
      addConnection(idx, idx_im1, vDy);
      addConnection(idx, idx_ip1, vDy);
      addConnection(idx, idx_jm1, vDx);
      addConnection(idx, idx_jp1, vDx);
    }
  }

  return derivedGeom;
}

// --------------------------------------------------------------------------------------

std::shared_ptr<Diffusion::DerivedGeom> Diffusion::calculateDerivedGeom(
    const oops::GeometryData & geom) {
  oops::Log::trace() << "Diffusion::calculateDerivedGeom start" << std::endl;
  util::Timer timer("oops::Diffusion", "calculateDerivedGeom");

  std::unique_ptr<Diffusion::DerivedGeom> derivedGeom;

  // create the derived geometry needed by diffusion, the method used depends on the atlas
  // functionspace type.
  const atlas::FunctionSpace & fs = geom.functionSpace();
  if (fs.type() == "NodeColumns") {
    derivedGeom = calculateDerivedGeom_NodeColumns(fs, geom);
  } else if (fs.type() == "StructuredColumns") {
    derivedGeom = calculateDerivedGeom_StructuredColumns(fs);
  } else {
    util::abor1_cpp("ERROR in oops::Diffusion(), atlas function type \""
                    +fs.type()+"\" is not supported");
  }

  oops::Log::trace() << "Diffusion::Diffusion end" << std::endl;
  return derivedGeom;
}

// --------------------------------------------------------------------------------------

Diffusion::Diffusion(const GeometryData & geometryData,
                     const std::shared_ptr<DerivedGeom> & derivedGeom)
  : geom_(geometryData), derivedGeom_(derivedGeom)
{
  oops::Log::trace() << "Diffusion::Diffusion start" << std::endl;
  util::Timer timer("oops::Diffusion", "Diffusion");

  if (derivedGeom_ == nullptr) {
    derivedGeom_ = Diffusion::calculateDerivedGeom(geom_);
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
  const std::vector<DerivedGeom::EdgeGeom> & edgeGeom = derivedGeom_->edgeGeom;
  parameters.haloExchange();

  // sanity check
  if (!parameters.has(HZ_SCALES) && !parameters.has(VT_SCALES)) {
    util::abor1_cpp("Diffusion::setScales() neither 'hzScales' nor 'vtScales' was given.",
      __FILE__, __LINE__);
  }

  // Horizontal diffusion parameters (in units of m)
  //-------------------------------------------------------------------------------------
  if (parameters.has(HZ_SCALES)) {
    // calculate the actual min number of hz iterations, and the diffusion coefficients (khdt_)
    double minItr = 0;
    const auto v_hzScales = atlas::array::make_view<double, 2>(parameters[HZ_SCALES]);
    khdtLevels_ = v_hzScales.shape(1);
    khdt_.resize(edgeGeom.size(), std::vector<double>(khdtLevels_));

    for (size_t e = 0; e < khdt_.size(); e++) {
      for (int level = 0; level < khdtLevels_; level++) {
        if (v_hzScales(edgeGeom[e].nodeA, level) == 0.0 ||
            v_hzScales(edgeGeom[e].nodeB, level) == 0.0) {
          // one of the nodes is masked out, dont diffuse with it
          khdt_[e][level] = 0.0;
        } else {
          // calculate diffusion coefficient (not taking into account the number of
          // iterations.. yet)
          const double s = (v_hzScales(edgeGeom[e].nodeA, level) +
                            v_hzScales(edgeGeom[e].nodeB, level)) / 2.0;
          khdt_[e][level] = s * s;

          // calculate the minimum number of iterations needed to be computationally stable
          // on this PE
          const double el = edgeGeom[e].edgeLength;
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
  const std::vector<DerivedGeom::EdgeGeom> & edgeGeom = derivedGeom_->edgeGeom;

  ASSERT(field.shape(0) == fs.size());
  ASSERT(khdtLevels_ == 1 || field.shape(1) <= khdtLevels_);

  const auto inv_area = atlas::array::make_view<double, 1>(derivedGeom_->inv_area);
  auto fieldVal = atlas::array::make_view<double, 2>(field);
  std::vector<double> flux(edgeGeom.size()*field.shape(1), 0.0);
  int fluxIdx;

  for (atlas::idx_t itr = 0; itr < niterHz_/2; itr++) {
    field.haloExchange();

    // NOTE: for efficiency with Atlas, the following loops are ordered here so
    // that we loop over the edges first, then the levels.

    // calculate the diffusive flux at each edge.
    fluxIdx = 0;
    for (size_t e = 0; e < edgeGeom.size(); e++) {
      const auto nodeA = edgeGeom[e].nodeA;
      const auto nodeB = edgeGeom[e].nodeB;
      const auto edgeAR = edgeGeom[e].aspectRatio;

      for (atlas::idx_t level = 0; level < field.shape(1); level++) {
        // khdtLevels is the level to pull khdt from, this code will work with
        // either a 3D khdt field or a 2D khdt field.
        const atlas::idx_t khdtLevel = std::min(level, khdtLevels_-1);
        const double dv = fieldVal(nodeA, level) - fieldVal(nodeB, level);
        flux[fluxIdx++] = (edgeAR * dv) * khdt_[e][khdtLevel];
      }
    }

    // time-step the diffusion terms.
    fluxIdx = 0;
      for (size_t e = 0; e < edgeGeom.size(); e++) {
      const auto nodeA = edgeGeom[e].nodeA;
      const auto nodeB = edgeGeom[e].nodeB;

      for (atlas::idx_t level = 0; level < field.shape(1); level++) {
        fieldVal(nodeA, level) -= inv_area(nodeA) * flux[fluxIdx];
        fieldVal(nodeB, level) += inv_area(nodeB) * flux[fluxIdx++];
      }
    }

    field.set_dirty(true);
  }
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplyHzAD(atlas::Field & field) const {
  const atlas::FunctionSpace & fs = geom_.functionSpace();
  const std::vector<DerivedGeom::EdgeGeom> & edgeGeom = derivedGeom_->edgeGeom;

  ASSERT(field.shape(0) == fs.size());
  ASSERT(khdtLevels_ == 1 || field.shape(1) <= khdtLevels_);

  const auto & inv_area = atlas::array::make_view<double, 1>(derivedGeom_->inv_area);
  const auto & ghost = atlas::array::make_view<int, 1>(fs.ghost());
  auto fieldVal = atlas::array::make_view<double, 2>(field);

  std::vector<double> flux(edgeGeom.size()*field.shape(1), 0.0);
  int fluxIdx;

  // init halo to zero
  for (atlas::idx_t i = 0; i < fs.size(); i++) {
    if (ghost(i)) {
      for (atlas::idx_t level = 0; level < field.shape(1); level++) {
        fieldVal(i, level) = 0;
      }
    }
  }

  for (int itr = 0; itr < niterHz_/2; itr++) {
    // NOTE: for efficiency with Atlas, the loops are ordered here so that we
    // loop over the edges first, then the levels.

    // adjoint time-step the diffusion terms
    fluxIdx = 0;
    for (size_t e = 0; e < edgeGeom.size(); e++) {
      const auto nodeA = edgeGeom[e].nodeA;
      const auto nodeB = edgeGeom[e].nodeB;
      for (int level = 0; level < field.shape(1); level++) {
        flux[fluxIdx++] = inv_area(nodeB) * fieldVal(nodeB, level)
                         -inv_area(nodeA) * fieldVal(nodeA, level);
      }
    }

    // Adjoint calculate diffusive flux at each edge.
    fluxIdx = 0;
    for (size_t e = 0; e < edgeGeom.size(); e++) {
      const auto nodeA = edgeGeom[e].nodeA;
      const auto nodeB = edgeGeom[e].nodeB;
      const auto edgeAR = edgeGeom[e].aspectRatio;

      for (int level = 0; level < field.shape(1); level++) {
        // khdtLevels is the level to pull khdt from, this code will work with
        // either a 3D khdt field or a 2D khdt field.
        const int khdtLevel = std::min(level, khdtLevels_-1);
        const auto dv = edgeAR * khdt_[e][khdtLevel] * flux[fluxIdx++];
        fieldVal(nodeA, level) += dv;
        fieldVal(nodeB, level) -= dv;
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
