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
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/generic/AtlasInterpolator.h"
#include "oops/generic/LocalInterpolatorBase.h"
#include "oops/generic/UnstructuredInterpolator.h"
#include "oops/runs/Test.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
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

// A direct test of the oops::AtlasInterpolator or oops::UnstructuredInterpolator, as just the
// interpolation matrix, without the communication of data to/from the interpolating processors.
// This is in contrast to how the interpolator would be used in practice, perhaps in GetValues or
// a grid change or a resolution change.
void testLocalInterpolator(
    const eckit::Configuration & config,
    const bool testSourcePointMask) {
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
  // can be selected via the config argument. Since we're not communicating across processors,
  // the target points must be drawn from a lon-lat region fully within the source grid. This
  // will be checked below.
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

  // Generate random points for the target grid
  util::UniformDistribution<double> rand_lons(nb_target, lon_bounds[0], lon_bounds[1], seed);
  util::UniformDistribution<double> rand_lats(nb_target, lat_bounds[0], lat_bounds[1], seed);
  const std::vector<double> target_lons = rand_lons.data();
  const std::vector<double> target_lats = rand_lats.data();

  std::unique_ptr<oops::LocalInterpolatorBase> interpolator;
  const eckit::LocalConfiguration interp_config(config, "local interpolator");
  const std::string interp_type = interp_config.getString("local interpolator type");
  if (interp_type == "atlas interpolator") {
    interpolator.reset(new oops::AtlasInterpolator(interp_config, geom, target_lats, target_lons));
  } else if (interp_type == "oops unstructured grid interpolator") {
    interpolator.reset(
        new oops::UnstructuredInterpolator(interp_config, geom, target_lats, target_lons));
  } else {
    throw eckit::BadParameter("Unknown local interpolator type: " + interp_type);
  }

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

  // If testing with masks, two additional tasks:
  // 1. add mask metadata to source_field
  // 2. add mask field to Geometry.fields()
  if (testSourcePointMask) {
    source_field.metadata().set("mask", "test_mask");

    // A simple mask for generic testing: mask = 0 in southern hemisphere; 1 in northern hemisphere
    atlas::Field mask = source_fs.createField<double>(name("test_mask") | levels(1));
    auto mask_view = make_view<double, 2>(mask);
    for (size_t jj = 0; jj < nb_source; ++jj) {
      const double lat = source_lonlat(jj, 1);
      mask_view(jj, 0) = (lat >= 0.0 ? 1.0 : 0.0);
    }
    // Hackily cast away the constness so we can shove a mask into the geometry FieldSet
    const_cast<atlas::FieldSet &>(geom.fieldSet()).add(mask);

    // We set the field values to a huge missingValue in the masked region, because this helps test
    // correctness at the boundaries of the mask. In detail: for stencils including masked and
    // unmasked source points, we still expect the interpolator to return a physical value, but if
    // there's a bug in the mask that allows a masked source point to (incorrectly) contribute to
    // the interpolation result, then the huge missing value will leak through and signal a bug.
    auto infield = make_view<double, 2>(source_field);
    for (size_t jj = 0; jj < nb_source; ++jj) {
      const double lat = source_lonlat(jj, 1);
      if (lat < 0.0) {
        for (size_t jlev = 0; jlev < nb_levels; ++jlev) {
          infield(jj, jlev) = util::missingValue<double>();
        }
      }
    }
  }

  const oops::Variables vars({varname});

  // Interpolation
  std::vector<double> target_vals;
  interpolator->apply(vars, source_fields, target_vals);

  // Test whether interpolation error falls within tolerance
  const double tolerance = config.getDouble("tolerance interpolation");
  if (!testSourcePointMask) {
    for (size_t jlev = 0; jlev < nb_levels; ++jlev) {
      for (size_t jj = 0; jj < nb_target; ++jj) {
        EXPECT(oops::is_close_absolute(
                 target_vals[jj + nb_target * jlev],
                 smooth_function(target_lons[jj], target_lats[jj], jlev, nb_levels),
                 tolerance));
      }
    }
  } else {
    // Applying an interpolation mask can reduce the number of source points within the
    // interpolation stencil, which can increase the interpolation error. So we loosen the test's
    // tolerance. To make sure that we still detect any errors in the masking, we set the source
    // field to huge missingValue above wherever it should be masked out, this way any errors
    // in the mask will contribute huge values and will exceed the test tolerance.
    const double masked_tol = 100.0 * tolerance;
    for (size_t jlev = 0; jlev < nb_levels; ++jlev) {
      for (size_t jj = 0; jj < nb_target; ++jj) {
        if (target_vals[jj + nb_target * jlev] == util::missingValue<double>()) {
          // Interpolation should only return 'missing' when all inputs are masked.
          // In this test, this should only happen in southern hemisphere.
          EXPECT(target_lats[jj] < 0.0);
        } else {
          EXPECT(oops::is_close_absolute(
                   target_vals[jj + nb_target * jlev],
                   smooth_function(target_lons[jj], target_lats[jj], jlev, nb_levels),
                   masked_tol));
        }
      }
    }
  }

  // Now test the adjoint

  // Set up a random field on the target grid
  util::UniformDistribution<double> rand_vals(nb_target * nb_levels, 0.0, 1.0);
  const std::vector<double> target_vals_ad = rand_vals.data();

  // Set up a FieldSet for the adjoint result on source grid
  atlas::FieldSet source_fields_ad;
  atlas::Field source_field_ad = source_fs.createField<double>(name(varname) | levels(nb_levels));
  auto source_ad_view = make_view<double, 2>(source_field_ad);
  source_ad_view.assign(0.0);  // initialize before adjoint accumulation
  if (testSourcePointMask) {
    source_field_ad.metadata().set("mask", "test_mask");
  }
  source_fields_ad.add(source_field_ad);

  // Apply interpolator adjoint
  interpolator->applyAD(vars, source_fields_ad, target_vals_ad);

  // Check adjoint with a manual dot product
  double dot1 = 0.0;
  double dot2 = 0.0;
  if (!testSourcePointMask) {
    for (size_t jlev = 0; jlev < nb_levels; ++jlev) {
      for (size_t jj = 0; jj < nb_source; ++jj) {
        dot1 += source_view(jj, jlev) * source_ad_view(jj, jlev);
      }
      for (size_t jj = 0; jj < nb_target; ++jj) {
        dot2 += target_vals[jj + nb_target * jlev] * target_vals_ad[jj + nb_target * jlev];
      }
    }
  } else {
    for (size_t jlev = 0; jlev < nb_levels; ++jlev) {
      for (size_t jj = 0; jj < nb_source; ++jj) {
        const double lat = source_lonlat(jj, 1);
        if (lat >= 0.0) {
          dot1 += source_view(jj, jlev) * source_ad_view(jj, jlev);
        }
      }
      for (size_t jj = 0; jj < nb_target; ++jj) {
        if (target_vals[jj + nb_target * jlev] != util::missingValue<double>()) {
          dot2 += target_vals[jj + nb_target * jlev]
                  * target_vals_ad[jj + nb_target * jlev];
        }
      }
    }
  }

  const double toleranceAD = config.getDouble("tolerance AD", 1.0e-11);
  EXPECT(oops::is_close(dot1, dot2, toleranceAD));
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

  // Test the interpolator with/without using an interpolation mask
  testLocalInterpolator(config, false);
  if (config.getBool("test with simple source-point mask", false)) {
    testLocalInterpolator(config, true);
  }
}

class LocalInterpolator : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~LocalInterpolator() {}

 private:
  std::string testid() const override {return "test::LocalInterpolator";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    const eckit::LocalConfiguration conf(::test::TestEnvironment::config());

    for (const std::string & testCaseName : conf.keys()) {
      const eckit::LocalConfiguration testCaseConf(conf, testCaseName);
      ts.emplace_back(
        CASE("generic/LocalInterpolator/" + testCaseName, testCaseConf) {
          testInterpolatorFromConfig(testCaseConf);
        });
    }
  }

  void clear() const override {}
};

}  // namespace test
