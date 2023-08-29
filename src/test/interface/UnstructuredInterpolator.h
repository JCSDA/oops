/*
 * (C) Copyright 2022 UCAR
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

#include <boost/noncopyable.hpp>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/grid.h"
#include "atlas/option.h"
#include "atlas/util/Point.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredInterpolator.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Random.h"

#include "test/interface/GeometryFixture.h"
#include "test/TestEnvironment.h"

using atlas::array::make_view;
using atlas::option::halo;
using atlas::option::levels;
using atlas::option::name;

namespace test {

/// Smooth function for testing interpolation
double testfunc(const double lon, const double lat, const size_t level, const size_t nlevels) {
  const double deg2rad = M_PI / 180.0;
  const double zz = static_cast<double>(level + 1) / static_cast<double>(nlevels);
  return (sin(deg2rad * lon) * cos(deg2rad * lat) + zz) / 1.5;
}

/// Test OOPS unstructured-grid interpolator
template <typename MODEL>
void testInterpolator(const bool testSourcePointMask) {
  typedef oops::Geometry<MODEL> Geometry_;

  const eckit::Configuration &config = ::test::TestEnvironment::config();

  // Parse inputs
  const size_t num_target = static_cast<size_t>(config.getInt("number of target points"));
  const std::vector<double> lat_bounds = config.getDoubleVector("latitude bounds of target points");
  const std::vector<double> lon_bounds =
      config.getDoubleVector("longitude bounds of target points");
  ASSERT(lat_bounds.size() == 2);
  ASSERT(lat_bounds[1] > lat_bounds[0]);
  ASSERT(lon_bounds.size() == 2);
  ASSERT(lon_bounds[1] > lon_bounds[0]);

  const unsigned int default_seed = 452938;
  const unsigned int seed = static_cast<unsigned int>(config.getInt("seed", default_seed));

  const size_t nlev = static_cast<size_t>(config.getInt("number of levels"));
  ASSERT(nlev > 0);  // to ensure the rank-2 atlas::Fields expected by OOPS
  const std::string varname = config.getString("variable name");

  const std::unique_ptr<const Geometry_> geom(new Geometry_(GeometryFixture<MODEL>::getParameters(),
                                                      oops::mpi::world(), oops::mpi::myself()));

  // Generate random point cloud for the target grid.
  // Here we want to test the interpolation matrix, without the communication of data to/from
  // the interpolating processors, as would be done in a global interpolation operation.
  // We do this by partitioning the target grid so that each processor's list of target points
  // is contained within its owned grid patch.
  const size_t my_task = geom->getComm().rank();
  const util::UniformDistribution<double> lons(num_target, lon_bounds[0], lon_bounds[1], seed);
  const util::UniformDistribution<double> lats(num_target, lat_bounds[0], lat_bounds[1], seed);
  atlas::PointXY point;
  std::vector<atlas::PointXY> target_points;
  std::vector<double> target_lats;
  std::vector<double> target_lons;
  for (size_t jj = 0; jj < num_target; ++jj) {
    if (geom->closestTask(lats[jj], lons[jj]) == my_task) {
      point.assign(lons[jj], lats[jj]);
      target_points.push_back(point);
      target_lats.push_back(lats[jj]);
      target_lons.push_back(lons[jj]);
    }
  }
  atlas::functionspace::PointCloud target_fs(target_points);

  oops::UnstructuredInterpolator interpolator(config, *geom, target_lats, target_lons);

  // Test print method
  oops::Log::info() << "Interpolator created:\n" << interpolator << std::endl;

  // Initialize the source FieldSet with a smooth function
  std::vector<double> source_lats{};
  std::vector<double> source_lons{};
  geom->latlon(source_lats, source_lons, true);  // include halo
  const size_t num_source = source_lats.size();
  const atlas::FunctionSpace & source_fs = geom->functionSpace();
  atlas::Field source_field = source_fs->createField<double>(name(varname) | levels(nlev));
  source_field.metadata().set("interp_type", "default");
  auto source_view = make_view<double, 2>(source_field);
  for (size_t jj = 0; jj < num_source; ++jj) {
    for (size_t jlev = 0; jlev < nlev; ++jlev) {
      source_view(jj, jlev) = testfunc(source_lons[jj], source_lats[jj], jlev, nlev);
    }
  }
  atlas::FieldSet source_fields;
  source_fields.add(source_field);

  // If testing with masks, two additional tasks:
  // 1. add interp_source_point_mask metadata to source_field
  // 2. add mask field to Geometry.fields() (we could use one of the "native" masks for MODEL,
  //    but this would make it very hard to write a generic test).
  if (testSourcePointMask) {
    source_field.metadata().set("interp_source_point_mask", "testmask");

    // A simple mask for generic testing: mask = 0 in southern hemisphere; 1 in northern hemisphere
    atlas::Field mask = source_fs->createField<double>(name("testmask") | levels(1));
    auto mask_view = make_view<double, 2>(mask);
    for (size_t jj = 0; jj < num_source; ++jj) {
      mask_view(jj, 0) = (source_lats[jj] >= 0.0 ? 1.0 : 0.0);
    }
    // Hackily cast away the constness so we can shove a mask into the geometry
    const_cast<atlas::FieldSet &>(geom->geometry().fields()).add(mask);

    // We set the field values to a huge missingValue in the masked region, because this helps test
    // correctness at the boundaries of the mask. In detail: for stencils including masked and
    // unmasked source points, we still expect the interpolator to return a physical value, but if
    // there's a bug in the mask that allows a masked source point to (incorrectly) contribute to
    // the interpolation result, then the huge missing value will leak through and signal a bug.
    auto infield = make_view<double, 2>(source_field);
    for (size_t jj = 0; jj < num_source; ++jj) {
      if (source_lats[jj] < 0.0) {
        for (size_t jlev = 0; jlev < nlev; ++jlev) {
          infield(jj, jlev) = util::missingValue<double>();
        }
      }
    }
  }

  oops::Variables vars{};
  vars.push_back(varname);

  // Apply interpolation
  std::vector<double> target_vals;
  interpolator.apply(vars, source_fields, target_vals);

  // Get test tolerance
  const size_t my_num_target = target_lons.size();
  const double tolerance = config.getDouble("tolerance interpolation");
  if (!testSourcePointMask) {
    for (size_t jlev = 0; jlev < nlev; ++jlev) {
      for (size_t jj = 0; jj < my_num_target; ++jj) {
        EXPECT(oops::is_close_absolute(
                 target_vals[jj + my_num_target * jlev],
                 testfunc(target_lons[jj], target_lats[jj], jlev, nlev),
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
    for (size_t jlev = 0; jlev < nlev; ++jlev) {
      for (size_t jj = 0; jj < my_num_target; ++jj) {
        if (target_vals[jj + my_num_target * jlev] == util::missingValue<double>()) {
          // Interpolation should only return 'missing' when all inputs are masked.
          // In this test, this should only happen in southern hemisphere.
          EXPECT(target_lats[jj] < 0.0);
        } else {
          EXPECT(oops::is_close_absolute(
                   target_vals[jj + my_num_target * jlev],
                   testfunc(target_lons[jj], target_lats[jj], jlev, nlev),
                   masked_tol));
        }
      }
    }
  }

  // Now test the adjoint

  // Define random field on target grid
  const util::UniformDistribution<double> random_field(my_num_target * nlev, 0.0, 1.0);
  std::vector<double> target_vals_ad(my_num_target * nlev);
  for (size_t jj = 0; jj < target_vals_ad.size(); ++jj) {
    target_vals_ad[jj] = random_field[jj];
  }

  // Initialize FieldSet for the result on source grid. Set it to 0 before testing adjoint,
  // because this increments the output rather than setting it.
  atlas::FieldSet source_fields_ad;
  atlas::Field source_field_ad = source_fs.createField<double>(name(varname) | levels(nlev));
  auto source_ad_view = make_view<double, 2>(source_field_ad);
  for (size_t jlev = 0; jlev < nlev; ++jlev) {
    for (size_t jj = 0; jj < num_source; ++jj) {
      source_ad_view(jj, jlev) = 0.0;
    }
  }
  if (testSourcePointMask) {
    source_field_ad.metadata().set("interp_source_point_mask", "testmask");
  }
  source_fields_ad.add(source_field_ad);

  // Apply interpolator adjoint
  interpolator.applyAD(vars, source_fields_ad, target_vals_ad);

  // Check adjoint with a manual dot product
  double dot1 = 0.0;
  double dot2 = 0.0;
  if (!testSourcePointMask) {
    for (size_t jlev = 0; jlev < nlev; ++jlev) {
      for (size_t jj = 0; jj < num_source; ++jj) {
        dot1 += source_view(jj, jlev) * source_ad_view(jj, jlev);
      }
      for (size_t jj = 0; jj < my_num_target; ++jj) {
        dot2 += target_vals[jj + my_num_target * jlev] * target_vals_ad[jj + my_num_target * jlev];
      }
    }
  } else {
    for (size_t jlev = 0; jlev < nlev; ++jlev) {
      for (size_t jj = 0; jj < num_source; ++jj) {
        if (source_lats[jj] >= 0.0) {
          dot1 += source_view(jj, jlev) * source_ad_view(jj, jlev);
        }
      }
      for (size_t jj = 0; jj < my_num_target; ++jj) {
        if (target_vals[jj + my_num_target * jlev] != util::missingValue<double>()) {
          dot2 += target_vals[jj + my_num_target * jlev]
                  * target_vals_ad[jj + my_num_target * jlev];
        }
      }
    }
  }

  const double toleranceAD = config.getDouble("tolerance AD", 1.0e-11);
  EXPECT(oops::is_close(dot1, dot2, toleranceAD));
}

/// Oops interpolation interface test
template <typename MODEL>
class UnstructuredInterpolator : public oops::Test {
 public:
  UnstructuredInterpolator() {}
  virtual ~UnstructuredInterpolator() {}

 private:
  std::string testid() const override {return "test::UnstructuredInterpolator";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("generic/UnstructuredInterpolator/testInterpolator")
      { testInterpolator<MODEL>(false); });
    ts.emplace_back(CASE("generic/UnstructuredInterpolator/testSourcePointMask")
      { testInterpolator<MODEL>(true); });
  }

  void clear() const override {}
};

}  // namespace test
