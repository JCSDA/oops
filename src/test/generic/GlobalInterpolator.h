/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#pragma once

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/option.h"
#include "atlas/util/Point.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/generic/GlobalInterpolator.h"
#include "oops/runs/Test.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "test/TestEnvironment.h"

using atlas::array::make_view;
using atlas::option::levels;
using atlas::option::name;

namespace test {

double smooth_function(const double lon, const double lat,
                       const size_t level, const size_t nb_levels) {
  const double deg2rad = M_PI / 180.0;
  const double zz = static_cast<double>(level + 1) / static_cast<double>(nb_levels);
  return (sin(deg2rad * lon) * cos(deg2rad * lat) + zz) / 1.5;
}

// A test of the oops::GlobalInterpolator, going from a source grid specified in the yaml to a
// random target grid.
void testGlobalInterpolator(
    const eckit::Configuration & config) {
  // Set up source grid
  const eckit::LocalConfiguration source_config(config, "source grid");
  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace functionspace{};
  atlas::FieldSet fieldset{};

  const eckit::mpi::Comm & comm = oops::mpi::world();
  util::setupFunctionSpace(comm, source_config, grid, partitioner, mesh, functionspace, fieldset);
  const oops::GeometryData geom(functionspace, fieldset, true, comm);

  // Options for the target points: the number of target points and their lon-lat bounding box
  // can be selected via the config argument.
  const eckit::LocalConfiguration target_config(config, "target points");
  const unsigned int seed = 452938;
  const size_t nb_target = static_cast<size_t>(target_config.getInt("number", 100));
  const std::vector<double> lon_bounds =
    target_config.getDoubleVector("longitude bounds", std::vector<double>({{0., 360.}}));
  const std::vector<double> lat_bounds =
    target_config.getDoubleVector("latitude bounds", std::vector<double>({{-90., 90.}}));
  ASSERT(lon_bounds.size() == 2);
  ASSERT(lon_bounds[1] > lon_bounds[0]);
  ASSERT(lat_bounds.size() == 2);
  ASSERT(lat_bounds[1] > lat_bounds[0]);

  // Generate a random point cloud for the target grid. We manually partition the target-point
  // domain across MPI ranks in a very simple way: we split the longitude direction in nb_task
  // slices, each of which covers the full latitude interval.
  const size_t nb_tasks = comm.size();
  const size_t my_task = comm.rank();
  const double delta_lon = (lon_bounds[1] - lon_bounds[0]) / nb_tasks;
  const double my_min_lon = lon_bounds[0] + my_task * delta_lon;
  const double my_max_lon = my_min_lon + delta_lon;
  const util::UniformDistribution<double> rand_lons(nb_target, my_min_lon, my_max_lon, seed);
  const util::UniformDistribution<double> rand_lats(nb_target, lat_bounds[0], lat_bounds[1], seed);
  std::vector<atlas::PointXY> target_points(nb_target);
  for (size_t jj = 0; jj < nb_target; ++jj) {
    target_points[jj] = atlas::PointXY(rand_lons[jj], rand_lats[jj]);
  }
  const atlas::functionspace::PointCloud target_fs(target_points);

  const eckit::LocalConfiguration interp_config(config, "global interpolator");
  const oops::GlobalInterpolator interpolator(interp_config, geom, target_fs, comm);

  // Set up the source FieldSet, containing one smooth function
  const std::string varname = "test_variable";
  const size_t nb_levels = 3;
  const auto source_fs = geom.functionSpace();
  const auto source_lonlat = make_view<double, 2>(source_fs.lonlat());
  atlas::Field source_field = source_fs.createField<double>(name(varname) | levels(nb_levels));
  source_field.metadata().set("interp_type", "default");
  const size_t nb_source = source_lonlat.shape(0);
  auto source_view = make_view<double, 2>(source_field);
  for (size_t jj = 0; jj < nb_source; ++jj) {
    const double lon = source_lonlat(jj, 0);
    const double lat = source_lonlat(jj, 1);
    for (size_t jlev = 0; jlev < nb_levels; ++jlev) {
      source_view(jj, jlev) = smooth_function(lon, lat, jlev, nb_levels);
    }
  }
  atlas::FieldSet source_fields;
  source_fields.add(source_field);

  // Interpolate
  atlas::FieldSet target_fields;
  interpolator.apply(source_fields, target_fields);

  // Test whether interpolation error falls within tolerance
  const double tolerance = config.getDouble("tolerance interpolation");
  auto target_view = make_view<double, 2>(target_fields.field(varname));
  for (size_t jj = 0; jj < nb_target; ++jj) {
    for (size_t jlev = 0; jlev < nb_levels; ++jlev) {
      EXPECT(oops::is_close_absolute(
               target_view(jj, jlev),
               smooth_function(rand_lons[jj], rand_lats[jj], jlev, nb_levels),
               tolerance));
    }
  }

  // Now test the adjoint

  // Set up a random field on the target grid
  const oops::Variables vars({varname});
  const std::vector<size_t> varSizes({nb_levels});
  const atlas::FieldSet target_fields_ad = util::createRandomFieldSet(
      comm, target_fs, varSizes, vars.variables());

  // Apply interpolator adjoint
  atlas::FieldSet source_fields_ad;
  interpolator.applyAD(source_fields_ad, target_fields_ad);

  // Check adjoint
  const double dot_source = util::dotProductFieldSets(
      source_fields, source_fields_ad, vars.variables(), comm);
  const double dot_target = util::dotProductFieldSets(
      target_fields, target_fields_ad, vars.variables(), comm);
  const double toleranceAD = config.getDouble("tolerance AD", 1.0e-11);
  EXPECT(oops::is_close(dot_source, dot_target, toleranceAD));
}

// Test interpolation from an input yaml
// Take config by non-const copy so we can make local modifications as needed
void testInterpolatorFromConfig(eckit::LocalConfiguration config) {
  // We allow a special key to request converting a named grid to an unstructured grid
  if (config.has("source grid.build unstructured grid from named grid")) {
    const eckit::LocalConfiguration conf(config,
        "source grid.build unstructured grid from named grid");
    atlas::Grid named_grid(conf);
    std::vector<double> vectorLonLat;
    for (auto lonlat : named_grid.lonlat()) {
      vectorLonLat.push_back(lonlat[0]);
      vectorLonLat.push_back(lonlat[1]);
    }
    config.set("source grid.grid.type", "unstructured");
    config.set("source grid.grid.xy", vectorLonLat);
  }

  testGlobalInterpolator(config);
}

class GlobalInterpolator : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~GlobalInterpolator() {}

 private:
  std::string testid() const override {return "test::GlobalInterpolator";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());

    for (const std::string & testCaseName : conf.keys()) {
      const eckit::LocalConfiguration testCaseConf(conf, testCaseName);
      ts.emplace_back(
        CASE("generic/GlobalInterpolator/" + testCaseName, testCaseConf) {
          testInterpolatorFromConfig(testCaseConf);
        });
    }
  }

  void clear() const override {}
};

}  // namespace test
