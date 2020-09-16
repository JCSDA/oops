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
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

void testConstructor() {
  std::unique_ptr<oops::Variables> vars(new oops::Variables());
  EXPECT(vars.get());

  const std::vector<std::string> varnames{"bt", "emiss"};
  const std::vector<int> channels{1, 2, 3, 4};
  std::unique_ptr<oops::Variables> other(new oops::Variables(varnames, channels));
  EXPECT(other.get());

  vars.reset();
  EXPECT(!vars.get());
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
  }

  void clear() const override {}
};

}  // namespace test

#endif  // TEST_BASE_VARIABLES_H_
