/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#pragma once

#include <cmath>
#include <iostream>
#include <limits>
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
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/util/Point.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/base/Variables.h"
#include "oops/generic/GlobalTemplatedInterpolator.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
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

/// Test simple oops unstructured-mesh interpolator
template <typename MODEL>
void testInterpolator() {
  typedef oops::Geometry<MODEL> Geometry_;

  const eckit::Configuration &conf = ::test::TestEnvironment::config();

  // Parse inputs
  const size_t num_target = static_cast<size_t>(
      conf.getInt("number of target points per mpi rank"));
  const std::vector<double> lat_bounds = conf.getDoubleVector("latitude bounds of target points");
  const std::vector<double> lon_bounds = conf.getDoubleVector("longitude bounds of target points");
  ASSERT(lat_bounds.size() == 2);
  ASSERT(lat_bounds[1] > lat_bounds[0]);
  ASSERT(lon_bounds.size() == 2);
  ASSERT(lon_bounds[1] > lon_bounds[0]);

  const unsigned int default_seed = 452938;
  const unsigned int seed = static_cast<unsigned int>(conf.getInt("seed", default_seed));
  const size_t nlev = static_cast<size_t>(conf.getInt("number of levels", 3));
  const std::string varname = conf.getString("variable name", "test field");

  const std::unique_ptr<const Geometry_> geom(new Geometry_(
        GeometryFixture<MODEL>::getParameters(), oops::mpi::world(), oops::mpi::myself()));

  // Generate a random point cloud for the target grid. We need each proc to own a distinct set
  // of target points, so we partition the requested lat/lon region across the available MPI tasks:
  const size_t num_tasks = oops::mpi::world().size();
  const size_t my_task = oops::mpi::world().rank();
  const double delta_lon = (lon_bounds[1] - lon_bounds[0]) / num_tasks;
  const double my_min_lon = lon_bounds[0] + my_task * delta_lon;
  const double my_max_lon = my_min_lon + delta_lon;
  const util::UniformDistribution<double> lons(num_target, my_min_lon, my_max_lon, seed);
  const util::UniformDistribution<double> lats(num_target, lat_bounds[0], lat_bounds[1], seed);
  atlas::PointXY point;
  std::vector<atlas::PointXY> target_points;
  for (size_t jj = 0; jj < num_target; ++jj) {
    point.assign(lons[jj], lats[jj]);
    target_points.push_back(point);
  }
  atlas::functionspace::PointCloud target_fs(target_points);

  oops::GlobalInterpolator<MODEL> interpolator(conf, *geom, target_fs, geom->getComm());

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
  const size_t rank = source_field.rank();
  ASSERT(rank == 1 || rank == 2);
  if (rank == 1) {
    auto source_view = make_view<double, 1>(source_field);
    for (size_t jj = 0; jj < num_source; ++jj) {
      source_view(jj) = testfunc(source_lons[jj], source_lats[jj], 0, nlev);
    }
  } else {
    auto source_view = make_view<double, 2>(source_field);
    for (size_t jj = 0; jj < num_source; ++jj) {
      for (size_t jlev = 0; jlev < nlev; ++jlev) {
        source_view(jj, jlev) = testfunc(source_lons[jj], source_lats[jj], jlev, nlev);
      }
    }
  }
  atlas::FieldSet source_fields;
  source_fields.add(source_field);

  // Initialize FieldSet for the result on target grid
  atlas::FieldSet target_fields;
  const atlas::Field target_field = target_fs.createField<double>(name(varname) | levels(nlev));
  target_fields.add(target_field);

  // Interpolate
  interpolator.apply(source_fields, target_fields);

  // Test whether interpolation error falls within tolerance
  const double tolerance = conf.getDouble("tolerance interpolation");
  if (rank == 1) {
    auto target_view = make_view<double, 1>(target_field);
    for (size_t jj = 0; jj < num_target; ++jj) {
      EXPECT(oops::is_close_absolute(
               target_view(jj),
               testfunc(lons[jj], lats[jj], 0, nlev),
               tolerance));
    }
  } else {
    auto target_view = make_view<double, 2>(target_field);
    for (size_t jlev = 0; jlev < nlev; ++jlev) {
      for (size_t jj = 0; jj < num_target; ++jj) {
        EXPECT(oops::is_close_absolute(
                 target_view(jj, jlev),
                 testfunc(lons[jj], lats[jj], jlev, nlev),
                 tolerance));
      }
    }
  }

  // Now test the adjoint

  // Initialize a random field on the target grid
  const util::UniformDistribution<double> random_field(num_target * nlev, 0.0, 1.0);
  atlas::Field target_field_ad = target_fs.createField<double>(name(varname) | levels(nlev));
  target_field_ad.metadata().set("interp_type", "default");
  if (rank == 1) {
    auto target_ad_view = make_view<double, 1>(target_field_ad);
    for (size_t jj = 0; jj < num_target; ++jj) {
      target_ad_view(jj) = random_field[jj];
    }
  } else {
    auto target_ad_view = make_view<double, 2>(target_field_ad);
    for (size_t jj = 0; jj < num_target; ++jj) {
      for (size_t jlev = 0; jlev < nlev; ++jlev) {
        target_ad_view(jj, jlev) = random_field[jj + num_target * jlev];
      }
    }
  }
  atlas::FieldSet target_fields_ad;
  target_fields_ad.add(target_field_ad);

  // Initialize FieldSet for the result on source grid. Set it to 0 before testing adjoint,
  // because this increments the output rather than setting it.
  atlas::FieldSet source_fields_ad;
  atlas::Field source_field_ad = source_fs->createField<double>(name(varname) | levels(nlev));
  if (rank == 1) {
    auto source_ad_view = make_view<double, 1>(source_field_ad);
    for (size_t jj = 0; jj < num_source; ++jj) {
      source_ad_view(jj) = 0.0;
    }
  } else {
    auto source_ad_view = make_view<double, 2>(source_field_ad);
    for (size_t jlev = 0; jlev < nlev; ++jlev) {
      for (size_t jj = 0; jj < num_source; ++jj) {
        source_ad_view(jj, jlev) = 0.0;
      }
    }
  }
  source_fields_ad.add(source_field_ad);

  // Apply interpolator adjoint
  interpolator.applyAD(source_fields_ad, target_fields_ad);

  // Check adjoint with a manual dot product
  double dot1 = 0.0;
  double dot2 = 0.0;
  if (rank == 1) {
    const auto source_view = make_view<double, 1>(source_field);
    const auto source_ad_view = make_view<double, 1>(source_field_ad);
    const auto target_ad_view = make_view<double, 1>(target_field_ad);
    const auto target_view = make_view<double, 1>(target_field);
    for (size_t jj = 0; jj < num_source; ++jj) {
      dot1 += source_view(jj) * source_ad_view(jj);
    }
    for (size_t jj = 0; jj < num_target; ++jj) {
      dot2 += target_view(jj) * target_ad_view(jj);
    }
  } else {
    const auto source_view = make_view<double, 2>(source_field);
    const auto source_ad_view = make_view<double, 2>(source_field_ad);
    const auto target_ad_view = make_view<double, 2>(target_field_ad);
    const auto target_view = make_view<double, 2>(target_field);
    for (size_t jlev = 0; jlev < nlev; ++jlev) {
      for (size_t jj = 0; jj < num_source; ++jj) {
        dot1 += source_view(jj, jlev) * source_ad_view(jj, jlev);
      }
      for (size_t jj = 0; jj < num_target; ++jj) {
        dot2 += target_view(jj, jlev) * target_ad_view(jj, jlev);
      }
    }
  }

  // Collect sums across MPI ranks. This is needed because each ranks owns unrelated subsets of
  // the source and target domains; for the dot products to match, the domains must match too.
  double dot1_sum = 0.0;
  double dot2_sum = 0.0;
  geom->getComm().allReduce(dot1, dot1_sum, eckit::mpi::Operation::Code::SUM);
  geom->getComm().allReduce(dot2, dot2_sum, eckit::mpi::Operation::Code::SUM);

  const double toleranceAD = conf.getDouble("tolerance AD", 1.0e-11);
  EXPECT(oops::is_close(dot1_sum, dot2_sum, toleranceAD));
}

///  Oops interpolation interface test
template <typename MODEL>
class GlobalInterpolator : public oops::Test {
 public:
  GlobalInterpolator() {}
  virtual ~GlobalInterpolator() {}

 private:
  std::string testid() const override {return "test::GlobalInterpolator";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("generic/GlobalInterpolator/testInterpolator")
      { testInterpolator<MODEL>(); });
  }

  void clear() const override {}
};

}  // namespace test
