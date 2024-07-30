/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "atlas/field.h"
#include "atlas/grid.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/FieldSets.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// \brief Tests FieldSets read ctor
void testFieldSets() {
  const eckit::Configuration & config = TestEnvironment::config();

  // Setup atlas FunctionSpace
  const auto & commGeom = oops::mpi::world();
  const auto & commTime = oops::mpi::myself();

  const auto fspaceConfig(config.getSubConfiguration("functionspace"));

  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace functionspace{};
  atlas::FieldSet fieldset{};

  util::setupFunctionSpace(commGeom, fspaceConfig, grid, partitioner,
                           mesh, functionspace, fieldset);

  // Read times
  const util::DateTime date(config.getString("date"));
  const std::vector<util::DateTime> dates({date});

  // Create two random oops FieldSet3D
  std::vector<std::string> varnames(config.getStringVector("variables"));
  std::vector<std::size_t> varSizes(config.getUnsignedVector("variable sizes"));
  oops::Variables vars;
  for (size_t ivar = 0; ivar < varnames.size(); ivar++) {
    eckit::LocalConfiguration conf;
    conf.set("levels", varSizes[ivar]);
    vars.push_back(oops::Variable(varnames[ivar], conf));
  }

  oops::Log::info() << "before init" << std::endl;
  oops::FieldSet3D fset1(date, commGeom);
  oops::Log::info() << "before rand" << std::endl;
  fset1.randomInit(functionspace, vars);
  oops::Log::info() << "before init" << std::endl;
  oops::FieldSet3D fset2(date, commGeom);
  oops::Log::info() << "before rand" << std::endl;
  fset2.randomInit(functionspace, vars);

  // Write to file
  const auto ioConfig = config.getSubConfiguration("input output members 1");
  const auto members = ioConfig.getSubConfigurations("members");
  ASSERT(members.size() == 2);
  oops::Log::info() << "before write" << std::endl;
  fset1.write(members[0]);
  fset2.write(members[1]);

  // Read back in using FieldSets() read ctor and list of members
  oops::FieldSets fsets(functionspace, vars, dates, ioConfig, commTime);

  // Compare fieldsets
  ASSERT(dates.size() == 1);
  const size_t it = 0;
  fsets(it, 0) -= fset1;
  EXPECT_EQUAL(fsets(it, 0).norm(vars), 0.0);
  fsets(it, 1) -= fset2;
  EXPECT_EQUAL(fsets(it, 1).norm(vars), 0.0);

  // Read back in using FieldSets() read ctor and member template
  const auto ioConfig2 = config.getSubConfiguration("input output members 2");
  oops::FieldSets fsets2(functionspace, vars, dates, ioConfig2, commTime);

  // Compare fieldsets
  fsets2(it, 0) -= fset1;
  EXPECT_EQUAL(fsets2(it, 0).norm(vars), 0.0);
  fsets2(it, 1) -= fset2;
  EXPECT_EQUAL(fsets2(it, 1).norm(vars), 0.0);
}
// -----------------------------------------------------------------------------

class FieldSets : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~FieldSets() = default;

 private:
  std::string testid() const override {return "test::FieldSets";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("base/FieldSets/testFieldSets")
      { testFieldSets(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
