/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

// Unit tests for oops::Diffusion: exercise the real operator through its
// public API on a small atlas geometry. Tests always run in
// Mode::VerticalOnly, so the horizontal mesh is inert and the per-column
// kernel produced by the vertical scheme (explicit or implicit) is what
// we measure. vtScales is uniform across every column, so column 0 is
// representative of the whole field.

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/option.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/runtime/Main.h"
#include "eckit/testing/Test.h"

#include "oops/base/GeometryData.h"
#include "oops/generic/Diffusion.h"
#include "oops/mpi/mpi.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

namespace {

// -----------------------------------------------------------------------------
// Minimal geometry scaffold. A small regular Gaussian grid gives us a real
// NodeColumns/StructuredColumns with edges, which is what Diffusion's
// constructor needs; the vertical tests never exercise those edges.
// -----------------------------------------------------------------------------
struct TinyGeom {
  atlas::Grid grid;
  atlas::grid::Partitioner partitioner;
  atlas::Mesh mesh;
  atlas::FunctionSpace functionSpace;
  atlas::FieldSet fieldset;
  std::unique_ptr<oops::GeometryData> geom;
};

std::unique_ptr<TinyGeom> makeTinyGeom() {
  auto g = std::make_unique<TinyGeom>();
  eckit::LocalConfiguration config;
  eckit::LocalConfiguration gridCfg;
  gridCfg.set("type", "regular_gaussian");
  gridCfg.set("N", 4);
  config.set("grid", gridCfg);
  config.set("function space", "StructuredColumns");
  config.set("halo", 1);
  util::setupFunctionSpace(oops::mpi::world(), config,
                           g->grid, g->partitioner, g->mesh,
                           g->functionSpace, g->fieldset);
  g->geom.reset(new oops::GeometryData(g->functionSpace, g->fieldset,
                                       true, oops::mpi::world()));
  return g;
}

// -----------------------------------------------------------------------------
// FieldSet helpers. Each test builds a fresh field so there's no alias
// between inputs to independent diffusion calls.
// -----------------------------------------------------------------------------

atlas::FieldSet makeDirac(const atlas::FunctionSpace & fs, int nz, int dirac) {
  atlas::FieldSet out;
  atlas::Field f = fs.createField<double>(atlas::option::name("x") | atlas::option::levels(nz));
  auto v = atlas::array::make_view<double, 2>(f);
  v.assign(0.0);
  for (atlas::idx_t i = 0; i < f.shape(0); ++i) v(i, dirac) = 1.0;
  out.add(f);
  return out;
}

atlas::FieldSet makeVtScales(const atlas::FunctionSpace & fs, int nz, double L) {
  atlas::FieldSet out;
  atlas::Field f = fs.createField<double>(atlas::option::name("vtScales")
                                          | atlas::option::levels(nz));
  auto v = atlas::array::make_view<double, 2>(f);
  v.assign(L);
  out.add(f);
  return out;
}

atlas::FieldSet makeSinusoid(const atlas::FunctionSpace & fs, int nz, double phase) {
  atlas::FieldSet out;
  atlas::Field f = fs.createField<double>(atlas::option::name("x") | atlas::option::levels(nz));
  auto v = atlas::array::make_view<double, 2>(f);
  for (atlas::idx_t i = 0; i < f.shape(0); ++i)
    for (int k = 0; k < nz; ++k)
      v(i, k) = std::sin(0.13 * k + phase + 0.01 * i);
  out.add(f);
  return out;
}

std::vector<double> columnZero(const atlas::FieldSet & fset) {
  const auto v = atlas::array::make_view<double, 2>(fset.field(0));
  std::vector<double> out(v.shape(1));
  for (int k = 0; k < v.shape(1); ++k) out[k] = v(0, k);
  return out;
}

// Inner product over owned (non-ghost) nodes only. multiplySqrt{TL,AD} does
// halo exchanges internally, so ghost values are not guaranteed to match
// between TL and AD outputs; including them breaks <x,Ay> = <A^T x, y>.
double innerProduct(const atlas::FieldSet & a, const atlas::FieldSet & b,
                    const atlas::FunctionSpace & fs) {
  const auto va = atlas::array::make_view<double, 2>(a.field(0));
  const auto vb = atlas::array::make_view<double, 2>(b.field(0));
  const auto ghost = atlas::array::make_view<int, 1>(fs.ghost());
  double s = 0.0;
  for (atlas::idx_t i = 0; i < va.shape(0); ++i) {
    if (ghost(i)) continue;
    for (atlas::idx_t k = 0; k < va.shape(1); ++k)
      s += va(i, k) * vb(i, k);
  }
  return s;
}

int argMaxAbs(const std::vector<double> & v) {
  int best = 0;
  double bestVal = std::abs(v[0]);
  for (size_t k = 1; k < v.size(); ++k) {
    if (std::abs(v[k]) > bestVal) { bestVal = std::abs(v[k]); best = static_cast<int>(k); }
  }
  return best;
}

double halfWidth(const std::vector<double> & v, int peak) {
  const double peakVal = v[peak];
  if (peakVal <= 0.0) return 0.0;
  for (int k = 1; k < static_cast<int>(v.size()); ++k) {
    const int kUp   = std::min<int>(v.size() - 1, peak + k);
    const int kDown = std::max(0, peak - k);
    if (v[kUp] / peakVal < 0.5 || v[kDown] / peakVal < 0.5) return static_cast<double>(k);
  }
  return static_cast<double>(v.size());
}

// -----------------------------------------------------------------------------
// Explicit vertical scheme
// -----------------------------------------------------------------------------

CASE("explicit: dirac kernel is symmetric, peaked, non-negative, mass-conserving") {
  auto tg = makeTinyGeom();
  const int nz = 51;
  const int dirac = nz / 2;

  double prevPeak = std::numeric_limits<double>::infinity();
  for (double L : {1.0, 2.0, 5.0, 10.0}) {
    oops::Diffusion diffusion(*tg->geom);
    diffusion.setParameters(makeVtScales(tg->functionSpace, nz, L));

    auto state = makeDirac(tg->functionSpace, nz, dirac);
    diffusion.multiply(state, oops::Diffusion::Mode::VerticalOnly);

    const auto x = columnZero(state);
    EXPECT_EQUAL(argMaxAbs(x), dirac);
    for (int k = 1; k < dirac; ++k) {
      EXPECT(std::abs(x[dirac + k] - x[dirac - k])
             < 1e-10 * std::max(std::abs(x[dirac + k]), 1e-30));
    }
    for (double v : x) EXPECT(v >= -1e-14);

    double sum = 0.0;
    for (double v : x) sum += v;
    EXPECT(std::abs(sum - 1.0) < 1e-9);

    oops::Log::info() << "explicit L=" << L << "  peak=" << x[dirac]
                      << "  sum=" << sum << std::endl;
    EXPECT(x[dirac] < prevPeak + 1e-14);
    prevPeak = x[dirac];
  }
}

CASE("explicit: multiplySqrtTL and multiplySqrtAD are mutual adjoints") {
  auto tg = makeTinyGeom();
  const int nz = 50;
  for (double L : {1.0, 2.0, 5.0, 10.0}) {
    oops::Diffusion diffusion(*tg->geom);
    diffusion.setParameters(makeVtScales(tg->functionSpace, nz, L));

    const auto x = makeSinusoid(tg->functionSpace, nz, 0.7);
    const auto y = makeSinusoid(tg->functionSpace, nz, -0.3);

    atlas::FieldSet TLy; util::copyFieldSet(y, TLy);
    diffusion.multiplySqrtTL(TLy, oops::Diffusion::Mode::VerticalOnly);

    atlas::FieldSet ADx; util::copyFieldSet(x, ADx);
    diffusion.multiplySqrtAD(ADx, oops::Diffusion::Mode::VerticalOnly);

    const double lhs = innerProduct(x, TLy, tg->functionSpace);
    const double rhs = innerProduct(ADx, y, tg->functionSpace);
    oops::Log::info() << "explicit L=" << L
                      << "  <x,TL y>=" << lhs
                      << "  <AD x,y>=" << rhs
                      << "  rel_diff=" << std::abs(lhs - rhs) / std::max(std::abs(lhs), 1.0)
                      << std::endl;
    EXPECT(std::abs(lhs - rhs) < 1e-12 * std::max(std::abs(lhs), 1.0));
  }
}

// -----------------------------------------------------------------------------
// Implicit vertical scheme
// -----------------------------------------------------------------------------

CASE("implicit: multiplySqrtTL and multiplySqrtAD are mutual adjoints") {
  auto tg = makeTinyGeom();
  const int nz = 50;
  for (double L : {1.0, 5.0, 20.0, 50.0}) {
    for (int M : {2, 4}) {
      oops::Diffusion diffusion(*tg->geom);
      diffusion.setParameters(makeVtScales(tg->functionSpace, nz, L),
                              oops::Diffusion::VerticalMethod::Implicit, M);

      const auto x = makeSinusoid(tg->functionSpace, nz, 0.7);
      const auto y = makeSinusoid(tg->functionSpace, nz, -0.3);

      atlas::FieldSet TLy; util::copyFieldSet(y, TLy);
      diffusion.multiplySqrtTL(TLy, oops::Diffusion::Mode::VerticalOnly);

      atlas::FieldSet ADx; util::copyFieldSet(x, ADx);
      diffusion.multiplySqrtAD(ADx, oops::Diffusion::Mode::VerticalOnly);

      const double lhs = innerProduct(x, TLy, tg->functionSpace);
      const double rhs = innerProduct(ADx, y, tg->functionSpace);
      oops::Log::info() << "implicit L=" << L << " M=" << M
                        << "  rel_diff=" << std::abs(lhs - rhs) / std::max(std::abs(lhs), 1.0)
                        << std::endl;
      EXPECT(std::abs(lhs - rhs) < 1e-12 * std::max(std::abs(lhs), 1.0));
    }
  }
}

CASE("implicit: dirac kernel is symmetric, peaked, and half-width grows with L") {
  auto tg = makeTinyGeom();
  const int nz = 51;
  const int dirac = nz / 2;
  const int M = 2;

  double prevHalfWidth = 0.0;
  for (double L : {1.0, 2.0, 5.0, 10.0, 20.0}) {
    oops::Diffusion diffusion(*tg->geom);
    diffusion.setParameters(makeVtScales(tg->functionSpace, nz, L),
                            oops::Diffusion::VerticalMethod::Implicit, M);

    auto state = makeDirac(tg->functionSpace, nz, dirac);
    diffusion.multiply(state, oops::Diffusion::Mode::VerticalOnly);

    const auto x = columnZero(state);
    EXPECT_EQUAL(argMaxAbs(x), dirac);
    for (int k = 1; k < dirac; ++k) {
      EXPECT(std::abs(x[dirac + k] - x[dirac - k])
             < 1e-10 * std::max(std::abs(x[dirac + k]), 1e-30));
    }
    for (double v : x) EXPECT(v >= -1e-14);

    const double hw = halfWidth(x, dirac);
    oops::Log::info() << "implicit L=" << L << "  peak=" << x[dirac]
                      << "  half_width_levels=" << hw << std::endl;
    EXPECT(hw >= prevHalfWidth);
    prevHalfWidth = hw;
  }
}

CASE("implicit: user L is the Daley length of the output kernel (M=2 Matern-1)") {
  // For M=2 the kernel is (1 + r/L) exp(-r/L). Verify the shape produced by
  // the implicit solve matches this at r = L, 2L, 3L.
  auto tg = makeTinyGeom();
  const int nz = 101;
  const int dirac = nz / 2;
  const double L = 10.0;
  const int M = 2;

  oops::Diffusion diffusion(*tg->geom);
  diffusion.setParameters(makeVtScales(tg->functionSpace, nz, L),
                          oops::Diffusion::VerticalMethod::Implicit, M);

  auto state = makeDirac(tg->functionSpace, nz, dirac);
  diffusion.multiply(state, oops::Diffusion::Mode::VerticalOnly);
  const auto x = columnZero(state);
  const double peak = x[dirac];

  struct { int offsetInL; double expected; } probes[] = {
    { 1, 2.0 / std::exp(1.0) },
    { 2, 3.0 / std::exp(2.0) },
    { 3, 4.0 / std::exp(3.0) },
  };
  for (const auto & p : probes) {
    const int r = p.offsetInL * static_cast<int>(L);
    const double actual = x[dirac + r] / peak;
    const double tol = 0.01 * p.offsetInL + 0.005;
    oops::Log::info() << "M=2 L=" << L << " r=" << r
                      << " kernel=" << actual
                      << " analytic=" << p.expected
                      << " rel_err=" << std::abs(actual - p.expected) / p.expected
                      << std::endl;
    EXPECT(std::abs(actual - p.expected) < tol * p.expected);
  }
}

CASE("implicit: kernel tail narrows (approaches Gaussian) as M grows at fixed L") {
  auto tg = makeTinyGeom();
  const int nz = 101;
  const int dirac = nz / 2;
  const double L = 5.0;

  double prevTailRatio = 1.0;
  for (int M : {2, 4, 8, 16}) {
    oops::Diffusion diffusion(*tg->geom);
    diffusion.setParameters(makeVtScales(tg->functionSpace, nz, L),
                            oops::Diffusion::VerticalMethod::Implicit, M);

    auto state = makeDirac(tg->functionSpace, nz, dirac);
    diffusion.multiply(state, oops::Diffusion::Mode::VerticalOnly);
    const auto x = columnZero(state);

    const double peak     = x[dirac];
    const double atOneL   = x[dirac + static_cast<int>(L)];
    const double atThreeL = x[dirac + static_cast<int>(3 * L)];
    const double tailRatio = atThreeL / atOneL;

    oops::Log::info() << "L=" << L << " M=" << M
                      << " peak=" << peak
                      << " x(L)/peak=" << atOneL / peak
                      << " x(3L)/x(L)=" << tailRatio << std::endl;
    EXPECT(peak > 0.0);
    EXPECT(atOneL > 0.0);
    EXPECT(atThreeL >= 0.0);
    EXPECT(tailRatio < prevTailRatio + 1e-12);
    prevTailRatio = tailRatio;
  }
}

CASE("implicit: L = nz still diffuses across the full column at M=2") {
  auto tg = makeTinyGeom();
  const int nz = 51;
  const int dirac = nz / 2;
  const int M = 2;
  const double L = nz;

  oops::Diffusion diffusion(*tg->geom);
  diffusion.setParameters(makeVtScales(tg->functionSpace, nz, L),
                          oops::Diffusion::VerticalMethod::Implicit, M);

  auto state = makeDirac(tg->functionSpace, nz, dirac);
  diffusion.multiply(state, oops::Diffusion::Mode::VerticalOnly);
  const auto x = columnZero(state);

  const double peak     = x[dirac];
  const double atTop    = x[0];
  const double atBottom = x[nz - 1];
  const double boundaryRatio = std::min(atTop, atBottom) / peak;
  oops::Log::info() << "L=" << L << " M=" << M
                    << " peak=" << peak
                    << " boundary/peak=" << boundaryRatio << std::endl;
  EXPECT(boundaryRatio > 0.5);
}

CASE("implicit: very large L (>= nz) remains finite and non-negative") {
  auto tg = makeTinyGeom();
  const int nz = 50;
  const int dirac = nz / 2;
  const int M = 2;

  for (double L : {50.0, 100.0, 500.0}) {
    oops::Diffusion diffusion(*tg->geom);
    diffusion.setParameters(makeVtScales(tg->functionSpace, nz, L),
                            oops::Diffusion::VerticalMethod::Implicit, M);

    auto state = makeDirac(tg->functionSpace, nz, dirac);
    diffusion.multiply(state, oops::Diffusion::Mode::VerticalOnly);
    const auto x = columnZero(state);

    double maxVal = x[0];
    double minVal = x[0];
    for (double v : x) {
      EXPECT(std::isfinite(v));
      maxVal = std::max(maxVal, v);
      minVal = std::min(minVal, v);
    }
    EXPECT(minVal >= -1e-14);
    EXPECT(maxVal <= 1.0);

    const double flatness = (maxVal > 0.0) ? minVal / maxVal : 0.0;
    oops::Log::info() << "L=" << L << " min/max(kernel)=" << flatness
                      << " max=" << maxVal << std::endl;
    if (L >= 10.0 * nz) EXPECT(flatness > 0.99);
  }
}

// -----------------------------------------------------------------------------
// Horizontal diffusion on an unstructured (NodeColumns) mesh. When no model
// "area" field is supplied, Diffusion estimates each node's control-volume
// area from the mesh (median-dual). These tests check that estimate and the
// resulting kernel.
// -----------------------------------------------------------------------------

struct UnstructuredGeom {
  atlas::Grid grid;
  atlas::grid::Partitioner partitioner;
  atlas::Mesh mesh;
  atlas::FunctionSpace functionSpace;
  atlas::FieldSet fieldset;
  std::unique_ptr<oops::GeometryData> geom;
};

// Build a NodeColumns geometry from scattered locations (lonlat laid out
// [lon0, lat0, lon1, lat1, ...]).
std::unique_ptr<UnstructuredGeom> makeUnstructuredGeom(const std::vector<double> & lonlat) {
  auto g = std::make_unique<UnstructuredGeom>();
  eckit::LocalConfiguration cfg;
  cfg.set("function space", "NodeColumns");
  cfg.set("grid.type", "unstructured");
  cfg.set("grid.xy", lonlat);
  cfg.set("partitioner", "equal_regions");
  cfg.set("no point on last task", true);
  util::setupFunctionSpace(oops::mpi::world(), cfg, g->grid, g->partitioner,
                           g->mesh, g->functionSpace, g->fieldset);
  g->geom.reset(new oops::GeometryData(g->functionSpace, g->fieldset,
                                       true, oops::mpi::world()));
  return g;
}

// Hex-packed (triangular) lattice; interior nodes have degree 6.
std::vector<double> hexLattice(int nx, int ny, double d, double lon0, double lat0) {
  std::vector<double> p;
  for (int j = 0; j < ny; ++j) {
    const double lat = lat0 + j * d * std::sqrt(3.0) / 2.0;
    const double off = (j % 2) ? d / 2.0 : 0.0;
    for (int i = 0; i < nx; ++i) { p.push_back(lon0 + off + i * d); p.push_back(lat); }
  }
  return p;
}

// A central hub ringed to degree 8 plus concentric rings, giving the spread of
// node degrees an irregular correlated-R obs network produces.
std::vector<double> irregularHub() {
  std::vector<double> p;
  auto add = [&](double lon, double lat) { p.push_back(lon); p.push_back(lat); };
  auto ring = [&](int n, double r, double phase) {
    for (int k = 0; k < n; ++k) {
      const double a = 2.0 * M_PI * (k + phase) / n;
      add(r * std::cos(a), r * std::sin(a));
    }
  };
  add(0.0, 0.0);
  ring(8, 2.0, 0.00); ring(14, 4.2, 0.30); ring(20, 7.0, 0.15);
  return p;
}

// Independent per-node reference: the median-dual area (sum of (cell area)/(n
// vertices) over incident cells), incident-cell count, incident-edge count, and
// node xyz, all in the same 3D chordal metric the operator uses.
struct MeshRef {
  std::vector<double> dualArea, x, y, z;
  std::vector<int> cells, edges;
};

MeshRef analyzeMesh(atlas::Mesh & mesh) {
  if (!mesh.nodes().has_field("xyz")) atlas::mesh::actions::BuildXYZField()(mesh);
  try { atlas::mesh::actions::build_edges(mesh); } catch (...) {}

  const int n = mesh.nodes().size();
  MeshRef r;
  r.dualArea.assign(n, 0.0); r.cells.assign(n, 0); r.edges.assign(n, 0);
  r.x.resize(n); r.y.resize(n); r.z.resize(n);
  const auto xyz = atlas::array::make_view<double, 2>(mesh.nodes().field("xyz"));
  for (int i = 0; i < n; ++i) { r.x[i] = xyz(i, 0); r.y[i] = xyz(i, 1); r.z[i] = xyz(i, 2); }

  const auto & c2n = mesh.cells().node_connectivity();
  for (atlas::idx_t c = 0; c < mesh.cells().size(); ++c) {
    const atlas::idx_t ncols = c2n.cols(c);
    if (ncols < 3) continue;
    const atlas::idx_t n0 = c2n(c, 0);
    double cellArea = 0.0;
    for (atlas::idx_t t = 1; t + 1 < ncols; ++t) {
      const atlas::idx_t a = c2n(c, t), b = c2n(c, t + 1);
      const double ux = r.x[a]-r.x[n0], uy = r.y[a]-r.y[n0], uz = r.z[a]-r.z[n0];
      const double vx = r.x[b]-r.x[n0], vy = r.y[b]-r.y[n0], vz = r.z[b]-r.z[n0];
      const double cx = uy*vz-uz*vy, cy = uz*vx-ux*vz, cz = ux*vy-uy*vx;
      cellArea += 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
    }
    const double share = cellArea / static_cast<double>(ncols);
    for (atlas::idx_t v = 0; v < ncols; ++v) {
      r.dualArea[c2n(c, v)] += share;
      r.cells[c2n(c, v)]++;
    }
  }

  const auto & e2n = mesh.edges().node_connectivity();
  for (atlas::idx_t e = 0; e < mesh.edges().size(); ++e) {
    r.edges[e2n(e, 0)]++;
    r.edges[e2n(e, 1)]++;
  }
  return r;
}

// The per-node area the operator actually computed (no area field supplied).
std::vector<double> operatorArea(oops::GeometryData & geom) {
  oops::Diffusion diffusion(geom);
  const auto v = atlas::array::make_view<double, 1>(diffusion.inverseArea());
  std::vector<double> a(v.shape(0));
  for (atlas::idx_t i = 0; i < v.shape(0); ++i) a[i] = (v(i) > 0.0) ? 1.0 / v(i) : 0.0;
  return a;
}

CASE("horizontal: estimated node area matches the median-dual reference on an irregular mesh") {
  auto g = makeUnstructuredGeom(irregularHub());
  const auto area = operatorArea(*g->geom);
  auto ref = analyzeMesh(g->mesh);

  int checked = 0, minDeg = 1000, maxDeg = 0;
  for (size_t i = 0; i < ref.dualArea.size(); ++i) {
    if (ref.dualArea[i] <= 0.0) continue;
    // the operator's per-node area equals the independent median-dual reference
    // regardless of node degree (the area is summed per incident cell, not per edge)
    EXPECT(std::abs(area[i] - ref.dualArea[i]) < 1e-6 * ref.dualArea[i]);
    minDeg = std::min(minDeg, ref.cells[i]);
    maxDeg = std::max(maxDeg, ref.cells[i]);
    ++checked;
  }
  oops::Log::info() << "[hub] nodes checked=" << checked
                    << " degree range=" << minDeg << ".." << maxDeg << std::endl;
  EXPECT(checked > 0);
  EXPECT(minDeg < maxDeg);                   // a real spread of node degrees
  EXPECT(ref.cells[0] >= 7);                 // hub reaches high degree
  EXPECT_EQUAL(ref.edges[0], ref.cells[0]);  // and is interior (#edges == #cells)
}

CASE("horizontal: dirac kernel is stable and matches the Gaussian at a resolved scale") {
  auto g = makeUnstructuredGeom(hexLattice(21, 21, 2.0, -20.0, -20.0));
  auto ref = analyzeMesh(g->mesh);
  const auto fs = g->functionSpace;
  const int n = static_cast<int>(ref.dualArea.size());
  auto interior = [&](int i) { return ref.edges[i] == ref.cells[i] && ref.cells[i] == 6; };
  auto dist = [&](int a, int b) {
    const double dx = ref.x[a]-ref.x[b], dy = ref.y[a]-ref.y[b], dz = ref.z[a]-ref.z[b];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  };

  // mean interior edge length -> a well-resolved (4-edge) length scale
  const auto & e2n = g->mesh.edges().node_connectivity();
  double esum = 0.0; int ecnt = 0;
  for (atlas::idx_t e = 0; e < g->mesh.edges().size(); ++e) {
    const int a = e2n(e, 0), b = e2n(e, 1);
    if (ref.cells[a] == 6 && ref.cells[b] == 6) { esum += dist(a, b); ++ecnt; }
  }
  const double L = 4.0 * esum / ecnt;

  // dirac at the interior node nearest the patch centre
  const auto lonlat = atlas::array::make_view<double, 2>(g->mesh.nodes().lonlat());
  int p = -1; double best = 1e30;
  for (int i = 0; i < n; ++i) {
    if (!interior(i)) continue;
    const double d = lonlat(i, 0)*lonlat(i, 0) + lonlat(i, 1)*lonlat(i, 1);
    if (d < best) { best = d; p = i; }
  }
  EXPECT(p >= 0);

  oops::Diffusion diffusion(*g->geom);
  atlas::Field hz = fs.createField<double>(atlas::option::name("hzScales")
                                           | atlas::option::levels(1));
  atlas::array::make_view<double, 2>(hz).assign(L);
  atlas::FieldSet scales; scales.add(hz);
  diffusion.setParameters(scales);

  atlas::Field fld = fs.createField<double>(atlas::option::name("f") | atlas::option::levels(1));
  auto v = atlas::array::make_view<double, 2>(fld);
  v.assign(0.0); v(p, 0) = 1.0;
  atlas::FieldSet x; x.add(fld);
  diffusion.multiply(x, oops::Diffusion::Mode::HorizontalOnly);
  const auto vv = atlas::array::make_view<double, 2>(x.field("f"));

  const double peak = vv(p, 0);
  EXPECT(peak > 0.0);
  double maxVal = peak, minVal = peak; int q = -1; double qErr = 1e30;
  for (int i = 0; i < n; ++i) {
    if (!interior(i)) continue;
    maxVal = std::max(maxVal, vv(i, 0));
    minVal = std::min(minVal, vv(i, 0));
    const double e = std::abs(dist(p, i) - L);
    if (i != p && e < qErr) { qErr = e; q = i; }
  }
  EXPECT_EQUAL(maxVal, peak);       // peak stays at the dirac
  EXPECT(minVal >= -0.05 * peak);   // stable: only small explicit ringing, no blow-up

  const double corr = vv(q, 0) / peak;
  const double target = std::exp(-0.5 * dist(p, q)*dist(p, q) / (L * L));
  oops::Log::info() << "[dirac] r/L=" << dist(p, q)/L << " corr=" << corr
                    << " target=" << target << std::endl;
  EXPECT(std::abs(corr - target) < 0.03);
}

}  // namespace

int main(int argc, char **argv) {
  eckit::Main::initialise(argc, argv);
  return eckit::testing::run_tests(argc, argv);
}
