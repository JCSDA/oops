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


atlas::functionspace::NodeColumns fixture_fs_nodecolumns_halo_appended() {
  auto g = atlas::Grid("O8");
  auto mg = atlas::MeshGenerator("structured", atlas::util::Config("ghost_at_end", true));
  auto m = mg.generate(g);
  return atlas::functionspace::NodeColumns(m, atlas::option::halo(1) | atlas::option::levels(2));
}

atlas::functionspace::StructuredColumns fixture_fs_structuredcolumns() {
  auto g = atlas::Grid("O8");
  return atlas::functionspace::StructuredColumns(
            g, atlas::option::halo(1) | atlas::option::levels(2));
}

atlas::functionspace::CubedSphereNodeColumns fixture_fs_cubedspherenodecolumns() {
  const auto g = atlas::Grid("CS-LFR-12");
  const auto mg = atlas::MeshGenerator(g.meshgenerator() | atlas::option::halo(1));
  const auto m = mg.generate(g);
  return atlas::functionspace::CubedSphereNodeColumns(m, atlas::option::levels(2));
}

CASE("Test for_each_value with 2D NodeColumns fields including halos") {
  auto fs = fixture_fs_nodecolumns_halo_appended();

  auto fieldA = fs.createField<double>();
  auto fieldB = fs.createField<double>();
  auto fieldC = fs.createField<double>();

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);

  viewA.assign(1.);
  viewB.assign(2.);

  util::for_each_value(
    util::IndexRange::include_halo,
      [](const double& a, const double& b, double& c) { c = a + b; },
      fieldA, fieldB, fieldC);

  for (atlas::idx_t j = 0; j < viewC.shape(0); ++j) {
    for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
      EXPECT_EQUAL(viewC(j, i), 3.);
    }
  }
}

CASE("Test for_each_column with 2D NodeColumn fields including halos") {
  auto fs = fixture_fs_nodecolumns_halo_appended();

  auto fieldA = fs.createField<double>();
  auto fieldB = fs.createField<double>();
  auto fieldC = fs.createField<double>();

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);

  viewA.assign(1.);
  viewB.assign(2.);

  using View = atlas::array::LocalView<double, 1>;
  using ConstView = atlas::array::LocalView<const double, 1>;

  util::for_each_column(
    util::IndexRange::include_halo,
    [](ConstView aCol, ConstView bCol, View cCol) {
        for (atlas::idx_t i = 0; i < aCol.shape(0); ++i) {
          cCol(i) = aCol(i) + bCol(i);
        }
    },
    fieldA, fieldB, fieldC);

  for (atlas::idx_t j = 0; j < viewC.shape(0); ++j) {
    for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
      EXPECT_EQUAL(viewC(j, i), 3.);
    }
  }
}

CASE("Test for_each_value with 2D StructuredColumn fields excluding halos") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = fs.createField<double>();
  auto fieldB = fs.createField<double>();
  auto fieldC = fs.createField<double>();

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);

  viewA.assign(1.);
  viewB.assign(2.);
  viewC.assign(-1.);

  util::for_each_value(
    util::IndexRange::exclude_halo,
    [](const double& a, const double& b, double& c) { c = a + b; },
    fieldA, fieldB, fieldC);

  const auto mesh_ghost = atlas::array::make_view<int, 1>(fs.ghost());

  for (atlas::idx_t j = 0; j < viewC.shape(0); ++j) {
    if (mesh_ghost(j) > 0) {
      for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
        EXPECT_EQUAL(viewC(j, i), -1.);
      }
    } else {
      for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
        EXPECT_EQUAL(viewC(j, i), 3.);
      }
    }
  }
}

CASE("Test for_each_column with 2D StructuredColumn fields excluding halos") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = fs.createField<double>();
  auto fieldB = fs.createField<double>();
  auto fieldC = fs.createField<double>();

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);

  viewA.assign(1.);
  viewB.assign(2.);
  viewC.assign(-1.);

  using View = atlas::array::LocalView<double, 1>;
  using ConstView = atlas::array::LocalView<const double, 1>;

  util::for_each_column(
    util::IndexRange::exclude_halo,
    [](ConstView aCol, ConstView bCol, View cCol) {
        for (atlas::idx_t i = 0; i < aCol.shape(0); ++i) {
          cCol(i) = aCol(i) + bCol(i);
        }
    },
    fieldA, fieldB, fieldC);

  const auto mesh_ghost = atlas::array::make_view<int, 1>(fs.ghost());

  for (atlas::idx_t j = 0; j < viewC.shape(0); ++j) {
    if (mesh_ghost(j) > 0) {
      for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
        EXPECT_EQUAL(viewC(j, i), -1.);
      }
    } else {
      for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
        EXPECT_EQUAL(viewC(j, i), 3.);
      }
    }
  }
}

CASE("Test for_each_value with 2D fields using serial backend") {
  auto fs = fixture_fs_structuredcolumns();

  auto fieldA = fs.createField<double>();
  auto fieldB = fs.createField<double>();
  auto fieldC = fs.createField<double>();

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);

  viewA.assign(1.);
  viewB.assign(2.);

  util::for_each_value(
    util::ExecutionPattern::serial,
    util::IndexRange::include_halo,
    [](const double& a, const double& b, double& c) { c = a + b; },
    fieldA, fieldB, fieldC);

  for (atlas::idx_t j = 0; j < viewC.shape(0); ++j) {
    for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
      EXPECT_EQUAL(viewC(j, i), 3.);
    }
  }
}

CASE("Test for_each_value with 2D CubedSphereNodeColumns fields excluding halos") {
  auto fs = fixture_fs_cubedspherenodecolumns();

  auto fieldA = fs.createField<double>();
  auto fieldB = fs.createField<double>();
  auto fieldC = fs.createField<double>();

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);

  viewA.assign(1.);
  viewB.assign(2.);
  viewC.assign(-1.);

  util::for_each_value(
    util::IndexRange::exclude_halo,
    [](const double& a, const double& b, double& c) { c = a + b; },
    fieldA, fieldB, fieldC);

  const auto mesh_ghost = atlas::array::make_view<int, 1>(fs.ghost());

  for (atlas::idx_t j = 0; j < viewC.shape(0); ++j) {
    if (mesh_ghost(j) > 0) {
      for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
        EXPECT_EQUAL(viewC(j, i), -1.);
      }
    } else {
      for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
        EXPECT_EQUAL(viewC(j, i), 3.);
      }
    }
  }
}

CASE("Test for_each_value with 2D NodeColumn constructed "
     "from CubedSphereNodeColumns fields excluding halos") {
  auto cs = fixture_fs_cubedspherenodecolumns();
  atlas::functionspace::NodeColumns fs(cs);

  auto fieldA = fs.createField<double>();
  auto fieldB = fs.createField<double>();
  auto fieldC = fs.createField<double>();

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);

  viewA.assign(1.);
  viewB.assign(2.);
  viewC.assign(-1.);

  util::for_each_value(
    util::IndexRange::exclude_halo,
    [](const double& a, const double& b, double& c) { c = a + b; },
    fieldA, fieldB, fieldC);

  const auto mesh_ghost = atlas::array::make_view<int, 1>(fs.ghost());

  for (atlas::idx_t j = 0; j < viewC.shape(0); ++j) {
    if (mesh_ghost(j) > 0) {
      for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
        EXPECT_EQUAL(viewC(j, i), -1.);
      }
    } else {
      for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
        EXPECT_EQUAL(viewC(j, i), 3.);
      }
    }
  }
}

CASE("Test for_each_value with 2D NodeColumn fields excluding halos") {
  auto fs = fixture_fs_nodecolumns_halo_appended();

  auto fieldA = fs.createField<double>();
  auto fieldB = fs.createField<double>();
  auto fieldC = fs.createField<double>();

  auto viewA = atlas::array::make_view<double, 2>(fieldA);
  auto viewB = atlas::array::make_view<double, 2>(fieldB);
  auto viewC = atlas::array::make_view<double, 2>(fieldC);

  viewA.assign(1.);
  viewB.assign(2.);
  viewC.assign(-1.);

  // TODO(adamsl): remove false branch when Atlas version >= 0.43.0
  if constexpr (util::details::has_halo_member_v<atlas::Field>) {
    util::for_each_value(
      util::IndexRange::exclude_halo,
      [](const double& a, const double& b, double& c) { c = a + b; },
      fieldA, fieldB, fieldC);

    const auto mesh_ghost = atlas::array::make_view<int, 1>(fs.ghost());

    for (atlas::idx_t j = 0; j < viewC.shape(0); ++j) {
      if (mesh_ghost(j) > 0) {
        for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
          EXPECT_EQUAL(viewC(j, i), -1.);
        }
      } else {
        for (atlas::idx_t i = 0; i < viewC.shape(1); ++i) {
          EXPECT_EQUAL(viewC(j, i), 3.);
        }
      }
    }
  } else {
    EXPECT_THROWS_AS(
      util::for_each_value(
        util::IndexRange::exclude_halo,
        [](const double& a, const double& b, double& c) { c = a + b; },
        fieldA, fieldB, fieldC),
      eckit::NotImplemented);
  }
}

int main(int argc, char **argv) {
  return eckit::testing::run_tests(argc, argv);
}
