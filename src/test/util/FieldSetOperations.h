/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_FIELDSETOPERATIONS_H_
#define TEST_UTIL_FIELDSETOPERATIONS_H_

#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FieldSetOperations.h"


namespace test {

CASE("util/FieldSetOperations/StructuredColumns") {
  // FunctionSpace
  eckit::LocalConfiguration gridConfig;
  gridConfig.set("type", "regular_gaussian");
  gridConfig.set("N", "10");
  atlas::Grid grid(gridConfig);
  atlas::grid::Partitioner partitioner("equal_regions");
  atlas::grid::Distribution distribution(grid, partitioner);
  atlas::functionspace::StructuredColumns fspace(grid, distribution, atlas::option::halo(1));

  // FieldSet
  atlas::FieldSet fset;

  // GeometryData
  oops::GeometryData geometryData(fspace, fset, true, oops::mpi::world());

  // Variables sizes
  std::vector<size_t> variableSizes = {2, 4};

  // Variables
  oops::Variables vars({"var1", "var2"});

  // Create random fields
  atlas::FieldSet fset1 = util::createRandomFieldSet(geometryData, variableSizes, vars);
  const double dp1 = util::dotProductFieldSets(fset1, fset1, vars, geometryData.comm());
  EXPECT(oops::is_close(dp1, 5640.50122292, 1.0e-12));
  atlas::FieldSet fset2 = util::createRandomFieldSet(geometryData, variableSizes, vars);
  const double dp2 = util::dotProductFieldSets(fset2, fset2, vars, geometryData.comm());
  EXPECT(oops::is_close(dp2, 5546.6627708858978, 1.0e-12));

  // Copy FieldSet
  atlas::FieldSet fset1copy = util::copyFieldSet(fset1);
  const double dp1copy = util::dotProductFieldSets(fset1copy, fset1copy, vars, geometryData.comm());
  EXPECT(oops::is_close(dp1copy, dp1, 1.0e-12));

  // Share Fields
  atlas::FieldSet fset1sh = util::shareFields(fset1);
  const double dp1sh = util::dotProductFieldSets(fset1sh, fset1sh, vars, geometryData.comm());
  EXPECT(oops::is_close(dp1sh, dp1, 1.0e-12));

  // Remove Fields from FieldSet
  atlas::FieldSet fset1woVar1 = util::copyFieldSet(fset1);
  util::removeFieldsFromFieldSet(fset1woVar1, oops::Variables({"var1"}));
  atlas::FieldSet fset1woVar2 = util::copyFieldSet(fset1);
  util::removeFieldsFromFieldSet(fset1woVar2, oops::Variables({"var2"}));
  const double dp1woVar1 = util::dotProductFieldSets(fset1woVar1, fset1woVar1,
    oops::Variables({"var2"}), geometryData.comm());
  const double dp1woVar2 = util::dotProductFieldSets(fset1woVar2, fset1woVar2,
    oops::Variables({"var1"}), geometryData.comm());
  EXPECT(oops::is_close(dp1woVar1+dp1woVar2, dp1, 1.0e-12));

  // Get grid UID
  std::string uid = util::getGridUid(geometryData.functionSpace());
  EXPECT(uid == "2734f1f878e2e047d290b3a578fc2927");
  EXPECT(uid == util::getGridUid(fset1));

  // Set data to zero
  atlas::FieldSet fset1zero = util::copyFieldSet(fset1);
  util::zeroFieldSet(fset1zero);
  const double dpzero = util::dotProductFieldSets(fset1zero, fset1zero, vars, geometryData.comm());
  EXPECT(dpzero == 0.0);

  // Add FieldSets
  atlas::FieldSet fset1add = util::copyFieldSet(fset1);
  util::addFieldSets(fset1add, fset1);
  const double dp1add = util::dotProductFieldSets(fset1add, fset1add, vars, geometryData.comm());
  EXPECT(oops::is_close(dp1add, 4.0*dp1, 1.0e-12));

  // Multiply FieldSet
  atlas::FieldSet fset1mul = util::copyFieldSet(fset1);
  util::multiplyFieldSet(fset1mul, 3.0);
  const double dp1mul = util::dotProductFieldSets(fset1mul, fset1mul, vars, geometryData.comm());
  EXPECT(oops::is_close(dp1mul, 9.0*dp1, 1.0e-12));

  // Multiply FieldSets
  atlas::FieldSet fset1sq = util::copyFieldSet(fset1);
  util::multiplyFieldSets(fset1sq, fset1);
  const double dp1sq = util::dotProductFieldSets(fset1sq, fset1sq, vars, geometryData.comm());
  EXPECT(oops::is_close(dp1sq, 17050.525704646639, 1.0e-12));

  // Divide FieldSets
  atlas::FieldSet fset1div = util::copyFieldSet(fset1);
  util::divideFieldSets(fset1div, fset2);
  const double dp1div = util::dotProductFieldSets(fset1div, fset1div, vars, geometryData.comm());
  EXPECT(oops::is_close(dp1div, 591417953.9186399, 1.0e-12));

  // FieldSet square-root
  atlas::FieldSet fset1sqrt = util::copyFieldSet(fset1sq);
  util::sqrtFieldSet(fset1sqrt);
  const double dp1sqrt = util::dotProductFieldSets(fset1sqrt, fset1sqrt, vars, geometryData.comm());
  EXPECT(oops::is_close(dp1sqrt, dp1, 1.0e-12));
}

class FieldSetOperations : public oops::Test {
 private:
  std::string testid() const override {return "test::FieldSetOperations";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_FIELDSETOPERATIONS_H_
