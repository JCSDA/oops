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

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/Configuration.h"
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

  {
    const std::vector<std::string> varnames{};
    const std::vector<int> channels{};
    oops::Variables other(TestEnvironment::config(), "empty variables");
    EXPECT(other.variables() == varnames);
    EXPECT(other.channels() == channels);
  }

  {
    // Fixture
    oops::Variables other(std::vector<std::string>({"var1", "var2", "var3"}));
    other.addMetaData("var2", "levels", 20);
    other.addMetaData("var3", "levels", 30);
    oops::Log::info() << "variables local config: " << other << std::endl;

    // Testing .variablesMetaData()
    int modelLevels(0);
    std::vector<std::string> confKeys(other.variablesMetaData().keys());
    std::vector<std::string> refKeys{"var2", "var3"};
    EXPECT(confKeys == refKeys);

    eckit::LocalConfiguration confOut(other.variablesMetaData());
    int i(0);
    for (const std::string s : other.variablesMetaData().keys()) {
      modelLevels = other.getLevels(s);
      EXPECT_EQUAL(modelLevels, (i + 2) * 10);
      ++i;
    }
  }

  vars.reset();
  EXPECT(!vars.get());

  // Test construction from an eckit::Configuration entry that does not exist
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

void testArithmeticOperators() {
  std::vector<std::string> varsStartStr{"var1", "var2", "var3"};
  std::vector<std::string> varsFinalStr{"var1", "var2"};
  std::vector<std::string> varsAddRemStr{"var3"};

  // Check on removing other Variables object
  oops::Variables varsStartRemoveVars(varsStartStr);
  oops::Variables varsFinalRemoveVars(varsFinalStr);
  oops::Variables varsRemove(varsAddRemStr);
  varsStartRemoveVars -= varsRemove;
  EXPECT(varsStartRemoveVars == varsFinalRemoveVars);

  // Check on removing other string
  oops::Variables varsStartRemoveStr(varsStartStr);
  oops::Variables varsFinalRemoveStr(varsFinalStr);
  varsStartRemoveStr -= "var3";
  EXPECT(varsStartRemoveStr == varsFinalRemoveStr);

  // Check on adding other Variables object
  oops::Variables varsStartAddVars(varsFinalStr);
  oops::Variables varsFinalAddVars(varsStartStr);
  oops::Variables varsAdd(varsAddRemStr);
  varsStartAddVars += varsAdd;
  EXPECT(varsStartAddVars == varsFinalAddVars);

  // Check we get exception if we subtract vars with channels
  oops::Variables varsStartChnnl(varsStartStr);
  oops::Variables varsWithChannels{varsAddRemStr, std::vector<int>{1}};
  EXPECT_THROWS_AS(varsStartChnnl -= varsWithChannels, eckit::NotImplemented);
}

// -----------------------------------------------------------------------------

void testMetaDataArithmeticOperators() {
  // Fixture
  oops::Variables vars(std::vector<std::string>{"var1", "var2", "var3"});
  for (auto & i : std::vector<int>{1, 2, 3}) {
    std::string var("var" + std::to_string(i));
    vars.addMetaData(var, "levels", i * 10);
  }
  oops::Variables varsCopy(vars);

  oops::Variables var1(std::vector<std::string>{"var1"});
  var1.addMetaData("var1", "levels", 10);

  oops::Log::info() << "vars lconf = " << vars << std::endl;
  vars -= var1;
  oops::Log::info() << "vars lconf minus var1 = " << vars << std::endl;

  // check -= string
  oops::Variables vars23(vars);
  oops::Variables vars123(varsCopy);
  vars123 -= std::string{"var1"};
  oops::Log::info() << "vars23 = " << vars23 << std::endl;
  oops::Log::info() << "vars123 = " << vars123 << std::endl;
  EXPECT(vars123 == vars23);

  vars += var1;
  oops::Log::info() << "vars adding back var1 = " << vars << std::endl;
  EXPECT(vars == varsCopy);
}

// -----------------------------------------------------------------------------
/// \brief tests Variables::operator== and operator!=
void testEquality() {
  oops::Variables abc({"a", "b", "c"});
  oops::Variables acb({"a", "c", "b"});
  oops::Variables ba({"b", "a"});

  EXPECT(abc == acb);
  EXPECT(!(abc != acb));
  EXPECT(!(ba == abc));
  EXPECT(acb != ba);
}

// -----------------------------------------------------------------------------
/// \brief tests Variables::operator== and operator!=
void testEqualityWithMetaData() {
  oops::Variables vars123(std::vector<std::string>{"var1", "var2", "var3"});
  for (auto & i : std::vector<int>{1, 2, 3}) {
    std::string var("var" + std::to_string(i));
    vars123.addMetaData(var, "levels", i * 10);
  }

  oops::Log::info() << "vars123 = " << vars123 << std::endl;

  oops::Variables vars213(std::vector<std::string>{"var2", "var1", "var3"});
  for (auto & i : std::vector<int>{2, 1, 3}) {
    std::string var("var" + std::to_string(i));
    vars213.addMetaData(var, "levels", i * 10);
  }

  oops::Log::info() << "vars213 = " << vars213 << std::endl;

  oops::Variables vars213SameMeta(std::vector<std::string>{"var2", "var1", "var3"});
  for (auto & i : std::vector<int>{2, 1, 3}) {
    std::string var("var" + std::to_string(i));
    vars213SameMeta.addMetaData(var, "levels", 30);
  }

  oops::Log::info() << "vars213 same meta = " << vars213SameMeta << std::endl;

  EXPECT(vars123 == vars213);
  EXPECT(!(vars123 != vars213));

  EXPECT(vars123 != vars213SameMeta);
  EXPECT(!(vars123 == vars213SameMeta));
}


// -----------------------------------------------------------------------------
/// \brief tests Variables::intersection (also uses operator=)
void testIntersection() {
  oops::Variables empty;
  oops::Variables acb({"a", "c", "b"});
  oops::Variables b({"b"});
  oops::Variables ba({"b", "a"});
  oops::Variables de({"d", "e"});
  oops::Variables de12({"d", "e"}, std::vector<int>{1, 2});
  oops::Variables ed21({"e", "d"}, std::vector<int>{2, 1});
  oops::Variables de1({"d", "e"}, std::vector<int>{1});
  oops::Variables d12({"d"}, std::vector<int>{1, 2});
  oops::Variables b1({"b"}, std::vector<int>{1});
  oops::Variables d1({"d"}, std::vector<int>{1});

  oops::Variables test = empty;
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
    ts.emplace_back(CASE("Variables/testArithmeticOperators")
      { testArithmeticOperators(); });
    ts.emplace_back(CASE("Variables/testMetaDataArithmeticOperators")
      { testMetaDataArithmeticOperators(); });
    ts.emplace_back(CASE("Variables/testEquality")
      { testEquality(); });
    ts.emplace_back(CASE("Variables/testEqualityWithMetaData")
      { testEqualityWithMetaData(); });
    ts.emplace_back(CASE("Variables/testIntersection")
      { testIntersection(); });
  }

  void clear() const override {}
};

}  // namespace test

#endif  // TEST_BASE_VARIABLES_H_
