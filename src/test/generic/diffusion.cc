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

}  // namespace

int main(int argc, char **argv) {
  eckit::Main::initialise(argc, argv);
  return eckit::testing::run_tests(argc, argv);
}
