/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#include "oops/util/for_each.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "eckit/testing/Test.h"

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

template<typename FunctionSpaceTy>
atlas::Field create_iota_field(FunctionSpaceTy& fs) {
  auto field = fs.template createField<double>();
  auto view = atlas::array::make_view<double, 2>(field);
  for (atlas::idx_t j = 0; j < view.shape(0); ++j) {
    for (atlas::idx_t i = 0; i < view.shape(1); ++i) {
      view(j, i) = static_cast<double>(i + (j * view.shape(1)));
    }
  }
  return field;
}

template<typename FunctionSpaceTy>
atlas::Field create_ones_field(FunctionSpaceTy& fs) {
  auto field = fs.template createField<double>();
  auto view = atlas::array::make_view<double, 2>(field);
  view.assign(1.0);
  return field;
}

template<typename FunctionSpaceTy>
atlas::Field create_zeros_field(FunctionSpaceTy& fs) {
  auto field = fs.template createField<double>();
  auto view = atlas::array::make_view<double, 2>(field);
  view.assign(0.0);
  return field;
}

void assert_all_values_equal(
  atlas::array::ArrayView<double, 2> expected_value,
  atlas::array::ArrayView<double, 2> view
) {
  for (atlas::idx_t j = 0; j < view.shape(0); ++j) {
    for (atlas::idx_t i = 0; i < view.shape(1); ++i) {
      EXPECT_EQUAL(view(j, i), expected_value(j, i));
    }
  }
}

void assert_owned_values_equal(
  atlas::array::ArrayView<double, 2> expected_value,
  atlas::array::ArrayView<double, 2> view,
  atlas::array::ArrayView<int, 1> ghost
) {
  for (atlas::idx_t j = 0; j < view.shape(0); ++j) {
    for (atlas::idx_t i = 0; i < view.shape(1); ++i) {
      if (ghost(j) == 0) {
        EXPECT_EQUAL(view(j, i), expected_value(j, i));
      }
    }
  }
}

void assert_ghost_values_equal(
  double expected_value,
  atlas::array::ArrayView<double, 2> view,
  atlas::array::ArrayView<int, 1> ghost
) {
  for (atlas::idx_t j = 0; j < view.shape(0); ++j) {
    for (atlas::idx_t i = 0; i < view.shape(1); ++i) {
      if (ghost(j) > 0) {
        EXPECT_EQUAL(view(j, i), expected_value);
      }
    }
  }
}

CASE("Test for_each_value with 2D NodeColumns fields including halos") {
  auto fs = fixture_fs_nodecolumns_halo_appended();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_value(
    util::IndexRange::include_halo,
      [](const double& a, const double& b, double& c) { c = a - b; },
      fieldA, fieldB, fieldC);

  assert_all_values_equal(viewTest, viewC);
}

CASE("Test for_each_column with 2D NodeColumn fields including halos") {
  auto fs = fixture_fs_nodecolumns_halo_appended();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  using View = atlas::array::LocalView<double, 1>;
  using ConstView = atlas::array::LocalView<const double, 1>;

  util::for_each_column(
    util::IndexRange::include_halo,
    [](ConstView aCol, ConstView bCol, View cCol) {
        for (atlas::idx_t i = 0; i < aCol.shape(0); ++i) {
          cCol(i) = aCol(i) - bCol(i);
        }
    },
    fieldA, fieldB, fieldC);

  assert_all_values_equal(viewTest, viewC);
}

CASE("Test for_each_value with 2D StructuredColumn fields excluding halos") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_value(
    util::IndexRange::exclude_halo,
    [](const double& a, const double& b, double& c) { c = a - b; },
    fieldA, fieldB, fieldC);

  const auto mesh_ghost = atlas::array::make_view<int, 1>(
    fieldC.functionspace().ghost());

  assert_owned_values_equal(viewTest, viewC, mesh_ghost);
  assert_ghost_values_equal(1., viewC, mesh_ghost);
}

CASE("Test for_each_column with 2D StructuredColumn fields excluding halos") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  using View = atlas::array::LocalView<double, 1>;
  using ConstView = atlas::array::LocalView<const double, 1>;

  util::for_each_column(
    util::IndexRange::exclude_halo,
    [](ConstView aCol, ConstView bCol, View cCol) {
        for (atlas::idx_t i = 0; i < aCol.shape(0); ++i) {
          cCol(i) = aCol(i) - bCol(i);
        }
    },
    fieldA, fieldB, fieldC);

  const auto mesh_ghost = atlas::array::make_view<int, 1>(
    fieldC.functionspace().ghost());

  assert_owned_values_equal(viewTest, viewC, mesh_ghost);
  assert_ghost_values_equal(1., viewC, mesh_ghost);
}

CASE("Test for_each_value with 2D fields using serial backend") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_value(
    util::ExecutionPattern::serial,
    util::IndexRange::include_halo,
    [](const double& a, const double& b, double& c) { c = a - b; },
    fieldA, fieldB, fieldC);

  assert_all_values_equal(viewTest, viewC);
}

CASE("Test for_each_value with 2D CubedSphereNodeColumns fields excluding halos") {
  auto fs = fixture_fs_cubedspherenodecolumns();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_value(
    util::IndexRange::exclude_halo,
    [](const double& a, const double& b, double& c) { c = a - b; },
    fieldA, fieldB, fieldC);

  const auto mesh_ghost = atlas::array::make_view<int, 1>(
    fieldC.functionspace().ghost());

  assert_owned_values_equal(viewTest, viewC, mesh_ghost);
  assert_ghost_values_equal(1., viewC, mesh_ghost);
}

CASE("Test for_each_value with 2D NodeColumn constructed "
     "from CubedSphereNodeColumns fields excluding halos") {
  const auto g = atlas::Grid("CS-LFR-12");
  const auto mg = atlas::MeshGenerator(g.meshgenerator() | atlas::option::halo(1));
  const auto m = mg.generate(g);
  auto cs = atlas::functionspace::CubedSphereNodeColumns(m, atlas::option::levels(2));
  atlas::functionspace::NodeColumns fs(cs);

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_value(
    util::IndexRange::exclude_halo,
    [](const double& a, const double& b, double& c) { c = a - b; },
    fieldA, fieldB, fieldC);

  const auto mesh_ghost = atlas::array::make_view<int, 1>(
    fieldC.functionspace().ghost());

  assert_owned_values_equal(viewTest, viewC, mesh_ghost);
  assert_ghost_values_equal(1., viewC, mesh_ghost);
}

CASE("Test for_each_value with 2D NodeColumn fields excluding halos") {
  auto fs = fixture_fs_nodecolumns_halo_appended();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_value(
    util::IndexRange::exclude_halo,
    [](const double& a, const double& b, double& c) { c = a - b; },
    fieldA, fieldB, fieldC);

  const auto mesh_ghost = atlas::array::make_view<int, 1>(
    fieldC.functionspace().ghost());

  assert_owned_values_equal(viewTest, viewC, mesh_ghost);
  assert_ghost_values_equal(1., viewC, mesh_ghost);
}

CASE("Test for_each_index with 2D StructuredColumn fields excluding halos") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_index(
    util::make_index_space_2d(util::IndexRange::exclude_halo, fieldC),
    [=](atlas::idx_t i, atlas::idx_t j) mutable {
      viewC(i, j) = viewA(i, j) - viewB(i, j);
    });

  const auto mesh_ghost = atlas::array::make_view<int, 1>(
    fieldC.functionspace().ghost());

  assert_owned_values_equal(viewTest, viewC, mesh_ghost);
  assert_ghost_values_equal(1., viewC, mesh_ghost);
}

CASE("Test for_each_index with 2D StructuredColumn fields excluding halos") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  using View = atlas::array::LocalView<double, 1>;
  using ConstView = atlas::array::LocalView<const double, 1>;

  util::for_each_index(
    util::ExecutionPattern::parallel,
    make_index_space_1d(util::IndexRange::exclude_halo, fieldC),
    [=](atlas::idx_t i) mutable {
        ConstView aCol = viewA.slice(i, atlas::array::Range::all());
        ConstView bCol = viewB.slice(i, atlas::array::Range::all());
        View cCol = viewC.slice(i, atlas::array::Range::all());
        for (atlas::idx_t j = 0; j < aCol.shape(0); ++j) {
          cCol(j) = aCol(j) - bCol(j);
        }
    });

  const auto mesh_ghost = atlas::array::make_view<int, 1>(
    fieldC.functionspace().ghost());

  assert_owned_values_equal(viewTest, viewC, mesh_ghost);
  assert_ghost_values_equal(1., viewC, mesh_ghost);
}

CASE("Test for_each_value with 2D fields using serial backend") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_index(
    util::ExecutionPattern::serial,
    util::make_index_space_2d(util::IndexRange::include_halo, fieldC),
    [=](atlas::idx_t i, atlas::idx_t j) mutable {
      viewC(i, j) = viewA(i, j) - viewB(i, j);
    });

  assert_all_values_equal(viewTest, viewC);
}

CASE("Test for_each_value with 3D fields using default backend") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = fs.createField<double>(atlas::option::name("vector") | atlas::option::variables(2));
  auto fieldB = fs.createField<double>(atlas::option::name("vector") | atlas::option::variables(2));
  auto fieldC = fs.createField<double>(atlas::option::name("vector") | atlas::option::variables(2));
  auto fieldTest = fs.createField<double>(
    atlas::option::name("vector") | atlas::option::variables(2));

  auto viewA = atlas::array::make_view<double, 3>(fieldA);
  auto viewB = atlas::array::make_view<double, 3>(fieldB);
  auto viewC = atlas::array::make_view<double, 3>(fieldC);
  auto viewTest = atlas::array::make_view<double, 3>(fieldTest);

  for (atlas::idx_t k = 0; k < viewA.shape(0); ++k) {
    for (atlas::idx_t j = 0; j < viewA.shape(1); ++j) {
      for (atlas::idx_t i = 0; i < viewA.shape(2); ++i) {
        double iota = i + (j * viewA.shape(2)) + (k * viewA.shape(1) * viewA.shape(2));
        viewA(k, j, i) = iota;
        viewB(k, j, i) = iota;
      }
    }
  }
  viewC.assign(1.0);
  viewTest.assign(0.0);

  util::for_each_index(
    util::make_index_space_3d(util::IndexRange::include_halo, fieldC),
    [=](atlas::idx_t i, atlas::idx_t j, atlas::idx_t k) mutable {
      viewC(i, j, k) = viewA(i, j, k) - viewB(i, j, k);
    });

  for (atlas::idx_t k = 0; k < viewA.shape(0); ++k) {
    for (atlas::idx_t j = 0; j < viewA.shape(1); ++j) {
      for (atlas::idx_t i = 0; i < viewA.shape(2); ++i) {
        EXPECT_EQUAL(viewC(i, j, k), viewTest(i, j, k));
      }
    }
  }
}

CASE("Test for_each_value with 3D fields using serial backend") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = fs.createField<double>(atlas::option::name("vector") | atlas::option::variables(2));
  auto fieldB = fs.createField<double>(atlas::option::name("vector") | atlas::option::variables(2));
  auto fieldC = fs.createField<double>(atlas::option::name("vector") | atlas::option::variables(2));
  auto fieldTest = fs.createField<double>(
    atlas::option::name("vector") | atlas::option::variables(2));

  auto viewA = atlas::array::make_view<double, 3>(fieldA);
  auto viewB = atlas::array::make_view<double, 3>(fieldB);
  auto viewC = atlas::array::make_view<double, 3>(fieldC);
  auto viewTest = atlas::array::make_view<double, 3>(fieldTest);

  for (atlas::idx_t k = 0; k < viewA.shape(0); ++k) {
    for (atlas::idx_t j = 0; j < viewA.shape(1); ++j) {
      for (atlas::idx_t i = 0; i < viewA.shape(2); ++i) {
        double iota = i + (j * viewA.shape(2)) + (k * viewA.shape(1) * viewA.shape(2));
        viewA(k, j, i) = iota;
        viewB(k, j, i) = iota;
      }
    }
  }
  viewC.assign(1.0);
  viewTest.assign(0.0);

  util::for_each_index(
    util::ExecutionPattern::serial,
    util::make_index_space_3d(util::IndexRange::include_halo, fieldC),
    [=](atlas::idx_t i, atlas::idx_t j, atlas::idx_t k) mutable {
      viewC(i, j, k) = viewA(i, j, k) - viewB(i, j, k);
    });

  for (atlas::idx_t k = 0; k < viewA.shape(0); ++k) {
    for (atlas::idx_t j = 0; j < viewA.shape(1); ++j) {
      for (atlas::idx_t i = 0; i < viewA.shape(2); ++i) {
        EXPECT_EQUAL(viewC(i, j, k), viewTest(i, j, k));
      }
    }
  }
}

CASE("Test for_each_value with 3D fields using parallel backend") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = fs.createField<double>(atlas::option::name("vector") | atlas::option::variables(2));
  auto fieldB = fs.createField<double>(atlas::option::name("vector") | atlas::option::variables(2));
  auto fieldC = fs.createField<double>(atlas::option::name("vector") | atlas::option::variables(2));
  auto fieldTest = fs.createField<double>(
    atlas::option::name("vector") | atlas::option::variables(2));

  auto viewA = atlas::array::make_view<double, 3>(fieldA);
  auto viewB = atlas::array::make_view<double, 3>(fieldB);
  auto viewC = atlas::array::make_view<double, 3>(fieldC);
  auto viewTest = atlas::array::make_view<double, 3>(fieldTest);

  for (atlas::idx_t k = 0; k < viewA.shape(0); ++k) {
    for (atlas::idx_t j = 0; j < viewA.shape(1); ++j) {
      for (atlas::idx_t i = 0; i < viewA.shape(2); ++i) {
        double iota = i + (j * viewA.shape(2)) + (k * viewA.shape(1) * viewA.shape(2));
        viewA(k, j, i) = iota;
        viewB(k, j, i) = iota;
      }
    }
  }
  viewC.assign(1.0);
  viewTest.assign(0.0);

  util::for_each_index(
    util::ExecutionPattern::parallel,
    util::make_index_space_3d(util::IndexRange::exclude_halo, fieldC),
    [=](atlas::idx_t i, atlas::idx_t j, atlas::idx_t k) mutable {
      viewC(i, j, k) = viewA(i, j, k) - viewB(i, j, k);
    });

  const auto mesh_ghost = atlas::array::make_view<int, 1>(
    fieldC.functionspace().ghost());

  for (atlas::idx_t k = 0; k < viewA.shape(0); ++k) {
    for (atlas::idx_t j = 0; j < viewA.shape(1); ++j) {
      for (atlas::idx_t i = 0; i < viewA.shape(2); ++i) {
        if (mesh_ghost(k) > 0) {
          EXPECT_EQUAL(1., viewC(k, j, i));
        } else {
          EXPECT_EQUAL(viewTest(k, j, i), viewC(k, j, i));
        }
      }
    }
  }
}

CASE("Test for_each_index with 2D Spectral fields with "
     "parallel execution pattern for a real-representation "
     "of the spectral coefficients") {
  auto fs = fixture_fs_spectral();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_index(
    util::ExecutionPattern::parallel,
    util::make_index_space_spectral_1d(fieldC),
    [=](atlas::idx_t i, atlas::idx_t n, atlas::idx_t m) mutable {
      for (int l = 0; l < viewC.shape(1); ++l) {
        viewC(i, l) = viewA(i, l) - viewB(i, l);
      }
    });

  assert_all_values_equal(viewTest, viewC);
}

CASE("Test for_each_index with 2D Spectral fields with "
     "serial execution pattern for a real-representation "
     "of the spectral coefficients") {
  auto fs = fixture_fs_spectral();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_index(
    util::ExecutionPattern::serial,
    util::make_index_space_spectral_1d(fieldC),
    [=](atlas::idx_t i, atlas::idx_t n, atlas::idx_t m) mutable {
      for (int l = 0; l < viewC.shape(1); ++l) {
        viewC(i, l) = viewA(i, l) - viewB(i, l);
      }
    });

  assert_all_values_equal(viewTest, viewC);
}

CASE("Test for_each_index with 2D Spectral fields with "
     "default execution pattern for a real-representation "
     "of the spectral coefficients") {
  auto fs = fixture_fs_spectral();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_index(
    util::make_index_space_spectral_1d(fieldC),
    [=](atlas::idx_t i, atlas::idx_t n, atlas::idx_t m) mutable {
      for (int l = 0; l < viewC.shape(1); ++l) {
        viewC(i, l) = viewA(i, l) - viewB(i, l);
      }
    });

  assert_all_values_equal(viewTest, viewC);
}

CASE("Test for_each_index with 2D Spectral fields with "
     "parallel execution pattern for a real-representation "
     "of the spectral coefficients") {
  auto fs = fixture_fs_spectral();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_index(
    util::ExecutionPattern::parallel,
    util::make_index_space_spectral_1d(fieldC),
    [=](atlas::idx_t i, atlas::idx_t n, atlas::idx_t m) mutable {
      for (int l = 0; l < viewC.shape(1); ++l) {
        viewC(i, l) = viewA(i, l) - viewB(i, l);
      }
    });

  assert_all_values_equal(viewTest, viewC);
}

CASE("Test for_each_index with 2D Spectral fields with "
     "serial execution pattern for a complex-representation "
     "of the spectral coefficients") {
  auto fs = fixture_fs_spectral();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_index(
    util::ExecutionPattern::serial,
    util::make_index_space_spectral_1d(fieldC),
    [=](atlas::idx_t i_R, atlas::idx_t i_Z, atlas::idx_t n, atlas::idx_t m) mutable {
      std::array<atlas::idx_t, 2> i{{i_R, i_Z}};
      for (int part = 0; part < 2; ++part) {
        for (int l = 0; l < viewC.shape(1); ++l) {
          viewC(i[part], l) = viewA(i[part], l) - viewB(i[part], l);
        }
      }
    });

  assert_all_values_equal(viewTest, viewC);
}

CASE("Test for_each_index with 2D Spectral fields with "
     "default execution pattern for a complex-representation "
     "of the spectral coefficients") {
  auto fs = fixture_fs_spectral();

  auto fieldA = create_iota_field(fs);
  auto fieldB = create_iota_field(fs);
  auto fieldC = create_ones_field(fs);
  auto fieldTest = create_zeros_field(fs);

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);
  auto viewTest = atlas::array::make_view<double, 2>(fieldTest);

  util::for_each_index(
    util::make_index_space_spectral_1d(fieldC),
    [=](atlas::idx_t i_R, atlas::idx_t i_Z, atlas::idx_t n, atlas::idx_t m) mutable {
      std::array<atlas::idx_t, 2> i{{i_R, i_Z}};
      for (int part = 0; part < 2; ++part) {
        for (int l = 0; l < viewC.shape(1); ++l) {
          viewC(i[part], l) = viewA(i[part], l) - viewB(i[part], l);
        }
      }
    });

  assert_all_values_equal(viewTest, viewC);
}

int main(int argc, char **argv) {
  return eckit::testing::run_tests(argc, argv);
}
