/*
 * (C) Copyright 2024- UCAR.
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

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsVariables.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

void testConstructor() {
  std::unique_ptr<oops::ObsVariables> vars(new oops::ObsVariables());
  EXPECT(vars.get());

  {
    const std::vector<std::string> varnames{"bt", "emiss"};
    const std::vector<int> channels{1, 2, 3, 4};
    std::unique_ptr<oops::ObsVariables> other(new oops::ObsVariables(varnames, channels));
    EXPECT(other.get());
    const std::vector<std::string> expectedObsVariables{"bt_1", "bt_2", "bt_3", "bt_4",
                                                     "emiss_1", "emiss_2", "emiss_3", "emiss_4"};
    EXPECT(other->variables() == expectedObsVariables);
    EXPECT(other->channels() == channels);
  }

  {
    const std::vector<std::string> varnames{"bt", "emiss"};
    const std::vector<int> channels{};
    std::unique_ptr<oops::ObsVariables> other(new oops::ObsVariables(varnames, channels));
    EXPECT(other.get());
    EXPECT(other->variables() == varnames);
    EXPECT(other->channels() == channels);
  }

  {
    const std::vector<std::string> varnames{};
    const std::vector<int> channels{};
    oops::ObsVariables other(TestEnvironment::config(), "empty variables");
    EXPECT(other.variables() == varnames);
    EXPECT(other.channels() == channels);
  }

  vars.reset();
  EXPECT(!vars.get());

  // Test construction from an eckit::Configuration entry that does not exist
  EXPECT_THROWS_AS(auto bad_vars = oops::ObsVariables(TestEnvironment::config(), "FOO"),
                   eckit::BadParameter);
}

// -----------------------------------------------------------------------------

void testCopyConstructor() {
  eckit::LocalConfiguration conf(TestEnvironment::config());
  std::unique_ptr<oops::ObsVariables> vars(new oops::ObsVariables(conf, "test variables"));
  EXPECT(vars.get());

  std::unique_ptr<oops::ObsVariables> other(new oops::ObsVariables(*vars));
  EXPECT(other.get());

  other.reset();
  EXPECT(!other.get());

  EXPECT(vars.get());
}

// -----------------------------------------------------------------------------

void testArithmeticOperators() {
  std::vector<std::string> varsStartStr{"var1", "var2", "var3"};
  std::vector<std::string> varsFinalStr{"var1", "var2"};
  std::vector<std::string> varsAddRemStr{"var3"};

  // Check on adding other ObsVariables object
  oops::ObsVariables varsStartAddVars(varsFinalStr);
  oops::ObsVariables varsFinalAddVars(varsStartStr);
  oops::ObsVariables varsAdd(varsAddRemStr);
  varsStartAddVars += varsAdd;
  EXPECT(varsStartAddVars == varsFinalAddVars);
}

// -----------------------------------------------------------------------------
/// \brief tests ObsVariables::operator== and operator!=
void testEquality() {
  oops::ObsVariables abc({"a", "b", "c"});
  oops::ObsVariables acb({"a", "c", "b"});
  oops::ObsVariables ba({"b", "a"});

  EXPECT(abc == acb);
  EXPECT(!(abc != acb));
  EXPECT(!(ba == abc));
  EXPECT(acb != ba);
}

// -----------------------------------------------------------------------------
/// \brief tests ObsVariables::intersection (also uses operator=)
void testIntersection() {
  oops::ObsVariables empty;
  oops::ObsVariables acb({"a", "c", "b"});
  oops::ObsVariables b({"b"});
  oops::ObsVariables ba({"b", "a"});
  oops::ObsVariables de({"d", "e"});
  oops::ObsVariables de12({"d", "e"}, std::vector<int>{1, 2});
  oops::ObsVariables ed21({"e", "d"}, std::vector<int>{2, 1});
  oops::ObsVariables de1({"d", "e"}, std::vector<int>{1});
  oops::ObsVariables d12({"d"}, std::vector<int>{1, 2});
  oops::ObsVariables b1({"b"}, std::vector<int>{1});
  oops::ObsVariables d1({"d"}, std::vector<int>{1});

  oops::ObsVariables test = empty;
  test.intersection(empty);
  EXPECT(test == empty);

  test = empty;
  test.intersection(acb);
  EXPECT(test == empty);

  test = acb;
  test.intersection(empty);
  EXPECT(test == empty);

  test = acb;
  test.intersection(b);
  EXPECT(test == b);

  test = b;
  test.intersection(acb);
  EXPECT(test == b);

  test = acb;
  test.intersection(ba);
  EXPECT(test == ba);

  test = ba;
  test.intersection(acb);
  EXPECT(test == ba);

  test = acb;
  test.intersection(de);
  EXPECT(test == empty);

  test = de;
  test.intersection(acb);
  EXPECT(test == empty);

  test = de12;
  test.intersection(b1);
  EXPECT(test == empty);

  test = de12;
  test.intersection(d1);
  EXPECT(test == d1);

  test = de12;
  test.intersection(ed21);
  EXPECT(test == de12);

  test = de12;
  test.intersection(de1);
  EXPECT(test == de1);

  test = de12;
  test.intersection(d12);
  EXPECT(test == d12);

  test = de12;
  EXPECT_THROWS(test.intersection(de));

  test = de12;
  test.intersection(empty);
  EXPECT(test == empty);
}

// -----------------------------------------------------------------------------

class ObsVariables : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~ObsVariables() {}

 private:
  std::string testid() const override {return "test::ObsVariables";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ObsVariables/testConstructor")
      { testConstructor(); });
    ts.emplace_back(CASE("ObsVariables/testCopyConstructor")
      { testCopyConstructor(); });
    ts.emplace_back(CASE("ObsVariables/testArithmeticOperators")
      { testArithmeticOperators(); });
    ts.emplace_back(CASE("ObsVariables/testEquality")
      { testEquality(); });
    ts.emplace_back(CASE("ObsVariables/testIntersection")
      { testIntersection(); });
  }

  void clear() const override {}
};

}  // namespace test
