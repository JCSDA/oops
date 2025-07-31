/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#include "oops/util/reduction.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "eckit/testing/Test.h"
#include "oops/util/for_each.h"
#include "oops/util/missingValues.h"

using util::missingValue;

auto fixture_fs_nodecolumns_halo_appended() {
  auto g = atlas::Grid("O8");
  auto mg = atlas::MeshGenerator("structured", atlas::util::Config("ghost_at_end", true));
  auto m = mg.generate(g);
  return atlas::functionspace::NodeColumns(
    m, atlas::option::halo(1) | atlas::option::levels(2));
}

auto fixture_fs_structuredcolumns() {
  auto g = atlas::Grid("O8");
  return atlas::functionspace::StructuredColumns(
    g, atlas::option::halo(1) | atlas::option::levels(2));
}

auto fixture_fs_cubedspherenodecolumns() {
  const auto g = atlas::Grid("CS-LFR-12");
  const auto mg = atlas::MeshGenerator(g.meshgenerator() | atlas::option::halo(1));
  const auto m = mg.generate(g);
  return atlas::functionspace::CubedSphereNodeColumns(m, atlas::option::levels(2));
}

auto fixture_fs_spectral() {
  atlas::idx_t truncation = 14;
  return atlas::functionspace::Spectral(truncation, atlas::option::levels(7));
}

template <typename FunctionSpaceTy>
atlas::Field create_zeros_field(const FunctionSpaceTy& fs) {
  auto field = fs.template createField<double>();
  auto view = atlas::array::make_view<double, 2>(field);
  view.assign(0.0);
  return field;
}

template <typename FunctionSpaceTy>
atlas::Field create_ones_field(const FunctionSpaceTy& fs) {
  auto field = fs.template createField<double>();
  auto view = atlas::array::make_view<double, 2>(field);
  view.assign(1.0);
  return field;
}

template <typename FunctionSpaceTy>
atlas::Field create_iota_field(const FunctionSpaceTy& fs) {
  auto field = fs.template createField<double>();
  auto view = atlas::array::make_view<double, 2>(field);
  for (atlas::idx_t j = 0; j < view.shape(0); ++j) {
    for (atlas::idx_t i = 0; i < view.shape(1); ++i) {
      view(j, i) = static_cast<double>(i + (j * view.shape(1)));
    }
  }
  return field;
}

double compute_expected_dot_product(const atlas::Field& field1, const atlas::Field& field2) {
  auto view1 = atlas::array::make_view<double, 2>(field1);
  auto view2 = atlas::array::make_view<double, 2>(field2);
  auto ghost = atlas::array::make_view<int, 1>(field1.functionspace().ghost());

  double expected = 0.0;
  for (atlas::idx_t j = 0; j < view1.shape(0); ++j) {
    if (ghost(j) == 0) {  // Only owned points
      for (atlas::idx_t i = 0; i < view1.shape(1); ++i) {
        expected += view1(j, i) * view2(j, i);
      }
    }
  }
  return expected;
}

double compute_expected_dot_product_missing_values(const atlas::Field& field1,
                                                   const atlas::Field& field2) {
  auto view1 = atlas::array::make_view<double, 2>(field1);
  auto view2 = atlas::array::make_view<double, 2>(field2);
  auto ghost = atlas::array::make_view<int, 1>(field1.functionspace().ghost());

  double expected = 0.0;
  for (atlas::idx_t j = 0; j < view1.shape(0); ++j) {
    if (ghost(j) == 0) {  // Only owned points
      for (atlas::idx_t i = 0; i < view1.shape(1); ++i) {
        if (view1(j, i) != missingValue<double>() &&
            view2(j, i) != missingValue<double>()) {
          expected += view1(j, i) * view2(j, i);
        }
      }
    }
  }
  return expected;
}

double compute_expected_dot_product_spectral(const atlas::Field& field1,
                                             const atlas::Field& field2) {
  auto view1 = atlas::array::make_view<double, 2>(field1);
  auto view2 = atlas::array::make_view<double, 2>(field2);

  const auto range = util::make_index_space_spectral_1d(field1);
  const int nb_zonal_wavenumbers{static_cast<int>(range.zonal_wavenumbers().size())};

  double expected = 0.0;
  atlas::idx_t j = 0;
  for (int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
    const int m = range.zonal_wavenumbers(jm);
    for (atlas::idx_t n = m; n <= range.truncation(); ++n) {
      if (m == 0) {
        for (atlas::idx_t i = 0; i < field1.shape(1); ++i) {
          expected += view1(j, i) * view2(j, i);
        }
        j += 2;
      } else {
        for (atlas::idx_t i = 0; i < field1.shape(1); ++i) {
          expected += 2.0 * view1(j, i) * view2(j, i);
          expected += 2.0 * view1(j + 1, i) * view2(j + 1, i);
        }
        j += 2;
      }
    }
  }
  return expected;
}

CASE("Test dot_product with zero field") {
  auto fs = fixture_fs_nodecolumns_halo_appended();

  auto fieldA = create_zeros_field(fs);
  auto fieldB = create_iota_field(fs);

  double expected = 0.0;
  double result = util::dot_product_on_task(fieldA, fieldB);

  EXPECT_EQUAL(result, expected);
}

CASE("Test dot_product with ones field") {
  auto fs = fixture_fs_nodecolumns_halo_appended();

  auto fieldA = create_ones_field(fs);
  auto fieldB = create_iota_field(fs);

  double expected = compute_expected_dot_product(fieldA, fieldB);
  double result = util::dot_product_on_task(fieldA, fieldB);

  EXPECT_EQUAL(result, expected);
}

CASE("Test dot_product with iota fields") {
  auto fs = fixture_fs_nodecolumns_halo_appended();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);

  double expected = compute_expected_dot_product(fieldA, fieldB);
  double result = util::dot_product_on_task(fieldA, fieldB);

  EXPECT_EQUAL(result, expected);
}

CASE("Test dot_product with fields containing missing values") {
  auto fs = fixture_fs_nodecolumns_halo_appended();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);

  // Set first few points to missing values
  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  for (atlas::idx_t j = 0; j < 3; ++j) {
    viewA(j, 0) = missingValue<double>();
    viewB(j, 1) = missingValue<double>();
  }

  double expected = compute_expected_dot_product_missing_values(fieldA, fieldB);
  double result = util::dot_product_on_task(fieldA, fieldB);

  EXPECT_EQUAL(result, expected);
}

CASE("Test dot_product with StructuredColumn fields") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);

  double expected = compute_expected_dot_product(fieldA, fieldB);
  double result = util::dot_product_on_task(fieldA, fieldB);

  EXPECT_EQUAL(result, expected);
}

CASE("Test dot_product with CubedSphereNodeColumns fields") {
  auto fs = fixture_fs_cubedspherenodecolumns();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);

  double expected = compute_expected_dot_product(fieldA, fieldB);
  double result = util::dot_product_on_task(fieldA, fieldB);

  EXPECT_EQUAL(result, expected);
}

CASE("Test dot_product using serial execution pattern") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);

  double expected = compute_expected_dot_product(fieldA, fieldB);
  double result = util::dot_product_on_task(util::ExecutionPattern::serial, fieldA, fieldB);

  EXPECT_EQUAL(result, expected);
}

CASE("Test dot_product using parallel execution pattern") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);

  double expected = compute_expected_dot_product(fieldA, fieldB);
  double result = util::dot_product_on_task(util::ExecutionPattern::parallel, fieldA, fieldB);

  EXPECT_EQUAL(result, expected);
}

CASE("Test dot_product with Spectral fields") {
  auto fs = fixture_fs_spectral();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);

  double expected = compute_expected_dot_product_spectral(fieldA, fieldB);
  double result = util::dot_product_on_task(fieldA, fieldB);

  EXPECT_EQUAL(result, expected);
}

CASE("Test dot_product with Spectral fields using serial execution pattern") {
  auto fs = fixture_fs_spectral();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);

  double expected = compute_expected_dot_product_spectral(fieldA, fieldB);
  double result = util::dot_product_on_task(util::ExecutionPattern::serial, fieldA, fieldB);

  EXPECT_EQUAL(result, expected);
}

CASE("Test dot_product with Spectral fields using parallel execution pattern") {
  auto fs = fixture_fs_spectral();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);

  double expected = compute_expected_dot_product_spectral(fieldA, fieldB);
  double result = util::dot_product_on_task(util::ExecutionPattern::parallel, fieldA, fieldB);

  EXPECT_EQUAL(result, expected);
}

int main(int argc, char **argv) {
  return eckit::testing::run_tests(argc, argv);
}
