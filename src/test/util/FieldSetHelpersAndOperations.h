/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"

namespace test {

CASE("util/FieldSetHelpersAndOperations/StructuredColumns") {
  // Communicator
  const eckit::mpi::Comm & comm = oops::mpi::world();

  // FunctionSpace
  eckit::LocalConfiguration gridConfig;
  gridConfig.set("type", "regular_gaussian");
  gridConfig.set("N", "10");
  atlas::Grid grid(gridConfig);
  atlas::grid::Partitioner partitioner("equal_regions");
  atlas::grid::Distribution distribution(grid, partitioner);
  atlas::functionspace::StructuredColumns fspace(grid, distribution, atlas::option::halo(1));

  // Variables
  std::vector<std::string> varnames({"var1", "var2"});
  eckit::LocalConfiguration metaml2, metaml4, variablesconf;
  metaml2.set("levels", 2);
  metaml4.set("levels", 4);
  variablesconf.set(varnames[0], metaml2);
  variablesconf.set(varnames[1], metaml4);
  oops::Variables vars(variablesconf, varnames);

  // Create random fields
  atlas::FieldSet fset1 = util::createRandomFieldSet(comm, fspace, vars);
  const double dp1 = util::dotProductFieldSets(fset1, fset1, vars.variables(), comm);
  double dp1_fields = 0.0;
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    dp1_fields += util::dotProductFields(fset1.field(vars[jvar]), fset1.field(vars[jvar]), comm);
  }
  EXPECT(oops::is_close(dp1, 5640.50122292, 1.0e-12));
  EXPECT(oops::is_close(dp1_fields, dp1, 1.0e-12));
  atlas::FieldSet fset2 = util::createRandomFieldSet(comm, fspace, vars);
  const double dp2 = util::dotProductFieldSets(fset2, fset2, vars.variables(), comm);
  EXPECT(oops::is_close(dp2, 5546.6627708858978, 1.0e-12));

  // Copy FieldSet
  atlas::FieldSet fset1copy = util::copyFieldSet(fset1);
  const double dp1copy = util::dotProductFieldSets(fset1copy, fset1copy, vars.variables(), comm);
  EXPECT(oops::is_close(dp1copy, dp1, 1.0e-12));

  // Share Fields
  atlas::FieldSet fset1sh = util::shareFields(fset1);
  const double dp1sh = util::dotProductFieldSets(fset1sh, fset1sh, vars.variables(), comm);
  EXPECT(oops::is_close(dp1sh, dp1, 1.0e-12));

  // Remove Fields from FieldSet
  atlas::FieldSet fset1woVar1 = util::copyFieldSet(fset1);
  util::removeFieldsFromFieldSet(fset1woVar1, {"var1"});
  atlas::FieldSet fset1woVar2 = util::copyFieldSet(fset1);
  util::removeFieldsFromFieldSet(fset1woVar2, {"var2"});
  const double dp1woVar1 = util::dotProductFieldSets(fset1woVar1, fset1woVar1, {"var2"}, comm);
  const double dp1woVar2 = util::dotProductFieldSets(fset1woVar2, fset1woVar2, {"var1"}, comm);
  EXPECT(oops::is_close(dp1woVar1+dp1woVar2, dp1, 1.0e-12));

  // Get grid UID
  std::string uid = util::getGridUid(fspace);
  EXPECT(uid == "2734f1f878e2e047d290b3a578fc2927");
  EXPECT(uid == util::getGridUid(fset1));

  // Set data to zero
  atlas::FieldSet fset1zero = util::copyFieldSet(fset1);
  util::zeroFieldSet(fset1zero);
  const double dpzero = util::dotProductFieldSets(fset1zero, fset1zero, vars.variables(), comm);
  EXPECT(dpzero == 0.0);

  // Add FieldSets
  atlas::FieldSet fset1add = util::copyFieldSet(fset1);
  util::addFieldSets(fset1add, fset1);
  const double dp1add = util::dotProductFieldSets(fset1add, fset1add, vars.variables(), comm);
  EXPECT(oops::is_close(dp1add, 4.0*dp1, 1.0e-12));

  // Subtract FieldSets
  atlas::FieldSet fset1sub = util::copyFieldSet(fset1);
  util::subtractFieldSets(fset1sub, fset1);
  const double dp1sub = util::dotProductFieldSets(fset1sub, fset1sub, vars.variables(), comm);
  EXPECT(dp1sub == 0.0);

  // Multiply FieldSet
  atlas::FieldSet fset1mul = util::copyFieldSet(fset1);
  util::multiplyFieldSet(fset1mul, 3.0);
  const double dp1mul = util::dotProductFieldSets(fset1mul, fset1mul, vars.variables(), comm);
  EXPECT(oops::is_close(dp1mul, 9.0*dp1, 1.0e-12));

  // Multiply FieldSets
  atlas::FieldSet fset1sq = util::copyFieldSet(fset1);
  util::multiplyFieldSets(fset1sq, fset1);
  const double dp1sq = util::dotProductFieldSets(fset1sq, fset1sq, vars.variables(), comm);
  EXPECT(oops::is_close(dp1sq, 17050.525704646639, 1.0e-12));

  // Divide FieldSets
  atlas::FieldSet fset1div = util::copyFieldSet(fset1);
  util::divideFieldSets(fset1div, fset2);
  const double dp1div = util::dotProductFieldSets(fset1div, fset1div, vars.variables(), comm);
  EXPECT(oops::is_close(dp1div, 591417953.9186399, 1.0e-12));

  // FieldSet square-root
  atlas::FieldSet fset1sqrt = util::copyFieldSet(fset1sq);
  util::sqrtFieldSet(fset1sqrt);
  const double dp1sqrt = util::dotProductFieldSets(fset1sqrt, fset1sqrt, vars.variables(), comm);
  EXPECT(oops::is_close(dp1sqrt, dp1, 1.0e-12));

  // Compare FieldSets
  EXPECT(util::compareFieldSets(fset1, fset1copy));
  EXPECT_NOT(util::compareFieldSets(fset1, fset2));

  // FieldSet norm
  const double norm1 = util::normFieldSet(fset1, vars.variables(), comm);
  EXPECT(oops::is_close(norm1, 69.620477228480709, 1.0e-12));

  // Create smooth Fieldset
  atlas::FieldSet smoothfset1 = util::createSmoothFieldSet(comm, fspace, vars);
  smoothfset1.haloExchange();
  atlas::FieldSet smoothfset2;

  // Write to file
  eckit::LocalConfiguration lconf;
  lconf.set("filepath", "./StructuredColumns_fset");
  util::writeFieldSet(comm, lconf, smoothfset1);

  // Read from file
  util::readFieldSet(comm, fspace, vars, lconf, smoothfset2);

  // Compare smooth FieldSets
  EXPECT(util::compareFieldSets(smoothfset1, smoothfset2));

  // Write to file (one file per task)
  lconf.set("one file per task", true);
  util::writeFieldSet(comm, lconf, smoothfset1);

  // Read from file (one file per task)
  util::readFieldSet(comm, fspace, vars, lconf, smoothfset2);

  // Compare smooth FieldSets
  EXPECT(util::compareFieldSets(smoothfset1, smoothfset2));
}

CASE("util/FieldSetHelpersAndOperations/NodeColumns") {
  // Communicator
  const eckit::mpi::Comm & comm = oops::mpi::world();

  // FunctionSpace
  eckit::LocalConfiguration gridConfig;
  gridConfig.set("name", "CS-LFR-15");
  atlas::Grid grid(gridConfig);
  atlas::Mesh mesh = atlas::MeshGenerator("cubedsphere_dual").generate(grid);
  atlas::functionspace::CubedSphereNodeColumns fspace(mesh);

  // Variables
  std::vector<std::string> varnames({"var1", "var2"});
  eckit::LocalConfiguration metaml2, metaml4, variablesconf;
  metaml2.set("levels", 2);
  metaml4.set("levels", 4);
  variablesconf.set(varnames[0], metaml2);
  variablesconf.set(varnames[1], metaml4);
  oops::Variables vars(variablesconf, varnames);

  // Create random fields
  atlas::FieldSet fset1 = util::createRandomFieldSet(comm, fspace, vars);
  const double dp1 = util::dotProductFieldSets(fset1, fset1, vars.variables(), comm);
  EXPECT(oops::is_close(dp1, 10135.544564309557, 1.0e-12));
  atlas::FieldSet fset2 = util::createRandomFieldSet(comm, fspace, vars);
  const double dp2 = util::dotProductFieldSets(fset2, fset2, vars.variables(), comm);
  EXPECT(oops::is_close(dp2, 10129.286186156089, 1.0e-12));

  // Get grid UID
  std::string uid = util::getGridUid(fspace);
  EXPECT(uid == "f277fb2f5d7f540eae4d0eba82f94eea");
  EXPECT(uid == util::getGridUid(fset1));

  // Write to file
  eckit::LocalConfiguration lconf;
  lconf.set("filepath", "./CubedSphere_fset");
  util::writeFieldSet(comm, lconf, fset1);

  // Read from file
  util::readFieldSet(comm, fspace, vars, lconf, fset2);

  // Compare smooth FieldSets
  EXPECT(util::compareFieldSets(fset1, fset2));

  // Write to file (one file per task)
  lconf.set("one file per task", true);
  util::writeFieldSet(comm, lconf, fset1);

  // Read from file (one file per task)
  util::readFieldSet(comm, fspace, vars, lconf, fset2);

  // Compare smooth FieldSets
  EXPECT(util::compareFieldSets(fset1, fset2));
}

CASE("util/FieldSetHelpersAndOperations/PointCloud") {
  // Communicator
  const eckit::mpi::Comm & comm = oops::mpi::world();

  // FunctionSpace
  eckit::LocalConfiguration gridConfig;
  gridConfig.set("type", "unstructured");
  const std::vector<double> lonlat({0, 90, 0,  -90, 180, -90, 180, 90});
  gridConfig.set("xy", lonlat);
  atlas::Grid grid(gridConfig);
  atlas::functionspace::PointCloud fspace(grid);

  // Variables
  std::vector<std::string> varnames({"var1", "var2"});
  eckit::LocalConfiguration metaml2, metaml4, variablesconf;
  metaml2.set("levels", 2);
  metaml4.set("levels", 4);
  variablesconf.set(varnames[0], metaml2);
  variablesconf.set(varnames[1], metaml4);
  oops::Variables vars(variablesconf, varnames);

  // Create random fields
  atlas::FieldSet fset1 = util::createRandomFieldSet(comm, fspace, vars);
  const double dp1 = util::dotProductFieldSets(fset1, fset1, vars.variables(), comm);
  EXPECT(oops::is_close(dp1, 16.260813650827739, 1.0e-12));
  atlas::FieldSet fset2 = util::createRandomFieldSet(comm, fspace, vars);
  const double dp2 = util::dotProductFieldSets(fset2, fset2, vars.variables(), comm);
  EXPECT(oops::is_close(dp2, 16.828513519408055, 1.0e-12));

  // Get grid UID
  std::string uid = util::getGridUid(fspace);
  EXPECT(uid == "2ac7a1ed6f655833eb61b9ececd5e5eb");
  EXPECT(uid == util::getGridUid(fset1));

  // Write to file (one file per task)
  eckit::LocalConfiguration lconf;
  lconf.set("filepath", "./PointCloud_fset");
  lconf.set("one file per task", true);
  util::writeFieldSet(comm, lconf, fset1);

  // Read from file (one file per task)
  util::readFieldSet(comm, fspace, vars, lconf, fset2);

  // Compare smooth FieldSets
  EXPECT(util::compareFieldSets(fset1, fset2));
}

class FieldSetHelpersAndOperations : public oops::Test {
 private:
  std::string testid() const override {return "test::FieldSetHelpersAndOperations";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test
