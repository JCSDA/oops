/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_BASE_VARIABLES_H_
#define TEST_BASE_VARIABLES_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Variables.h"
#include "oops/runs/Test.h"
#include "test/base/Fortran.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

void testConstructor() {
  std::unique_ptr<oops::Variables> vars(new oops::Variables());
  EXPECT(vars.get());

  {
    const std::vector<std::string> varnames{"bt", "emiss"};
    const std::vector<int> channels{1, 2, 3, 4};
    std::unique_ptr<oops::Variables> other(new oops::Variables(varnames, channels));
    EXPECT(other.get());
    const std::vector<std::string> expectedVariables{"bt_1", "bt_2", "bt_3", "bt_4",
                                                     "emiss_1", "emiss_2", "emiss_3", "emiss_4"};
    EXPECT(other->variables() == expectedVariables);
    EXPECT(other->channels() == channels);
  }

  {
    const std::vector<std::string> varnames{"bt", "emiss"};
    const std::vector<int> channels{};
    std::unique_ptr<oops::Variables> other(new oops::Variables(varnames, channels));
    EXPECT(other.get());
    EXPECT(other->variables() == varnames);
    EXPECT(other->channels() == channels);
  }

  vars.reset();
  EXPECT(!vars.get());

  // Test bad construction from eckit::Configuration
  EXPECT_THROWS_AS(auto bad_vars = oops::Variables(TestEnvironment::config(), "FOO"),
                   eckit::BadParameter);
}

// -----------------------------------------------------------------------------

void testCopyConstructor() {
  eckit::LocalConfiguration conf(TestEnvironment::config());
  std::unique_ptr<oops::Variables> vars(new oops::Variables(conf, "test variables"));
  EXPECT(vars.get());

  std::unique_ptr<oops::Variables> other(new oops::Variables(*vars));
  EXPECT(other.get());

  other.reset();
  EXPECT(!other.get());

  EXPECT(vars.get());
}

// -----------------------------------------------------------------------------

void testFortranInterface() {
  oops::Variables vars;

  test_vars_interface_f(TestEnvironment::config(), vars);

  // test content of variable list
  std::vector<std::string> vars_check = TestEnvironment::config().getStringVector("test variables");

  // The fortran routine tests the push_back_vector method, with the variables
  // from the config file, as well as the push of a single variable name.
  // So, if both were successful, vars should contain one extra item in
  // the variable list
  EXPECT(vars.size() ==  vars_check.size()+1);

  for (std::size_t jvar = 0; jvar < vars_check.size(); ++jvar) {
    EXPECT(vars[jvar] == vars_check[jvar]);
  }
}

// -----------------------------------------------------------------------------

class Variables : public oops::Test {
 public:
  Variables() {}
  virtual ~Variables() {}
 private:
  std::string testid() const override {return "test::Variables";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("Variables/testConstructor")
      { testConstructor(); });
    ts.emplace_back(CASE("Variables/testCopyConstructor")
      { testCopyConstructor(); });
    ts.emplace_back(CASE("Variables/testFortranInterface")
      { testFortranInterface(); });
  }

  void clear() const override {}
};

}  // namespace test

#endif  // TEST_BASE_VARIABLES_H_
