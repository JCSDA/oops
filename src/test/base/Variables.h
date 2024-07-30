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
    const std::vector<std::string> varnames{"air_temperature", "air_pressure"};
    std::unique_ptr<oops::Variables> other(new oops::Variables(varnames));
    EXPECT(other.get());
    EXPECT(other->variables() == varnames);
  }

  {
    const std::vector<std::string> varnames{};
    oops::Variables other(TestEnvironment::config(), "empty variables");
    EXPECT(other.variables() == varnames);
  }

  {
    // Fixture
    oops::Variable var1("var1");
    eckit::LocalConfiguration conf2;
    conf2.set("levels", 20);
    conf2.set("test_double", 2.0);
    conf2.set("test_string", "2");
    oops::Variable var2("var2", conf2);
    oops::Variable var3("var3",
                        oops::VariableMetaData(oops::VerticalStagger::INTERFACE,
                                               oops::ModelDataType::Int32),
                        30);
    oops::Variables other({var1, var2, var3});
    oops::Log::info() << "variables local config: " << other << std::endl;
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
  std::vector<std::string> vars_check_string
                                      = TestEnvironment::config().getStringVector("test variables");
  oops::Variables vars_check_variables(vars_check_string);

  // The fortran routine tests the push_back_vector method, with the variables
  // from the config file, as well as the push of a single variable name.
  // So, if both were successful, vars should contain one extra item in
  // the variable list
  EXPECT(vars.size() ==  vars_check_variables.size()+1);

  for (std::size_t jvar = 0; jvar < vars_check_variables.size(); ++jvar) {
    EXPECT(vars[jvar].name() == vars_check_variables[jvar].name());
  }
}

// -----------------------------------------------------------------------------

void testVariableConstructorAndEqualsComparison() {
  std::vector<std::string> varsStartStr{"var1", "var2", "var3"};
  oops::Variable var1NoLevels("var1");
  EXPECT(var1NoLevels.name() == "var1");
  EXPECT(var1NoLevels.getLevels() == -1);
  EXPECT(var1NoLevels.stagger() == oops::defaultVerticalStagger);
  EXPECT(var1NoLevels.dataType() == oops::defaultDataType);

  oops::Variable var1WithSetMetadata("var1",
                                     oops::VariableMetaData(oops::VerticalStagger::INTERFACE,
                                                            oops::ModelDataType::Int32));
  EXPECT(var1WithSetMetadata.getLevels() == -1);
  EXPECT(var1WithSetMetadata.stagger() == oops::VerticalStagger::INTERFACE);
  EXPECT(var1WithSetMetadata.dataType() == oops::ModelDataType::Int32);

  oops::Variable var1With10Levels("var1", oops::VariableMetaData(), 10);
  EXPECT(var1With10Levels.getLevels() == 10);
  EXPECT(var1With10Levels.stagger() == oops::defaultVerticalStagger);
  EXPECT(var1With10Levels.dataType() == oops::defaultDataType);

  oops::Variable var1With20Levels("var1", eckit::LocalConfiguration().set("levels", 20));
  EXPECT(var1With20Levels.getLevels() == 20);
  EXPECT(var1With20Levels.stagger() == oops::defaultVerticalStagger);
  EXPECT(var1With20Levels.dataType() == oops::defaultDataType);

  oops::Variable var2("var2");

  EXPECT(var1NoLevels.metaData() != var1WithSetMetadata.metaData());
  EXPECT(var1With10Levels.metaData() == var1With20Levels.metaData());
  EXPECT(var1NoLevels != var1WithSetMetadata);
  EXPECT(var1NoLevels == var1With10Levels);  // levels are ignored in this comparison
  EXPECT(var1With10Levels != var1With20Levels);  // levels are not ignored in this comparison
  EXPECT(var1NoLevels != var2);
}

// -------------------------------------------------------------------------------------------------

void testPushBack() {
  oops::Variables vars;
  oops::Variable var1("var1");
  oops::Variable var2("var2");
  oops::Variable var3("var3");
  oops::Variable var1WithSetMetadata("var1",
                                     oops::VariableMetaData(oops::VerticalStagger::INTERFACE,
                                                            oops::ModelDataType::Int32));
  oops::Variable var2With10Levels("var2", oops::VariableMetaData(), 10);

  vars.push_back(var1);
  vars.push_back(var2);
  vars.push_back(var3);
  vars.push_back(var1WithSetMetadata);

  EXPECT(vars.size() == 4);
  EXPECT(vars[0].name() == "var1");
  EXPECT(vars[1].name() == "var2");
  EXPECT(vars[2].name() == "var3");
  EXPECT(vars[3].name() == "var1");
  EXPECT(vars[1].getLevels() == -1);

  // Test that adding the same variable twice does not change the size
  vars.push_back(var1);
  EXPECT(vars.size() == 4);
  // Test that adding the "same" variable with levels set updates the levels
  vars.push_back(var2With10Levels);
  EXPECT(vars.size() == 4);
  EXPECT(vars[1].getLevels() == 10);
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
  varsStartRemoveStr -= oops::Variable{"var3"};
  EXPECT(varsStartRemoveStr == varsFinalRemoveStr);

  // Check on adding other Variables object
  oops::Variables varsStartAddVars(varsFinalStr);
  oops::Variables varsFinalAddVars(varsStartStr);
  oops::Variables varsAdd(varsAddRemStr);
  varsStartAddVars += varsAdd;
  EXPECT(varsStartAddVars == varsFinalAddVars);
}

// -----------------------------------------------------------------------------

void testMetaDataArithmeticOperators() {
  // Fixture
  oops::Variables vars;
  for (const auto & i : std::vector<int>{1, 2, 3}) {
    const std::string var("var" + std::to_string(i));
    eckit::LocalConfiguration conf;
    conf.set("levels", i * 10);
    vars.push_back(oops::Variable(var, conf));
  }

  oops::Variables varsCopy(vars);

  oops::Variables var1;
  eckit::LocalConfiguration conf;
  conf.set("levels", 10);
  var1.push_back(oops::Variable("var1", conf));

  oops::Log::info() << "vars lconf = " << vars << std::endl;
  vars -= var1;
  oops::Log::info() << "vars lconf minus var1 = " << vars << std::endl;

  // check -= string
  oops::Variables vars23(vars);
  oops::Variables vars123(varsCopy);

  eckit::LocalConfiguration conf1;
  conf1.set("levels", 10);
  vars123 -= oops::Variable{"var1", conf1};
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
  oops::Variables abc(std::vector<std::string>{"a", "b", "c"});
  oops::Variables acb(std::vector<std::string>{"a", "c", "b"});
  oops::Variables ba(std::vector<std::string>{"b", "a"});

  EXPECT(abc == acb);
  EXPECT(!(abc != acb));
  EXPECT(!(ba == abc));
  EXPECT(acb != ba);
}

// -----------------------------------------------------------------------------
/// \brief tests Variables::operator== and operator!=
void testEqualityWithMetaData() {
  oops::Variables vars123;
  for (const auto & i : std::vector<int>{1, 2, 3}) {
    const std::string var("var" + std::to_string(i));
    eckit::LocalConfiguration conf;
    conf.set("levels", i * 10);
    vars123.push_back(oops::Variable(var, conf));
  }

  oops::Log::info() << "vars123 = " << vars123 << std::endl;

  oops::Variables vars213;
  for (const auto & i : std::vector<int>{2, 1, 3}) {
    const std::string var("var" + std::to_string(i));
    eckit::LocalConfiguration conf;
    conf.set("levels", i * 10);
    vars213.push_back(oops::Variable(var, conf));
  }

  oops::Log::info() << "vars213 = " << vars213 << std::endl;

  oops::Variables vars213SameMeta;
  for (auto & i : std::vector<int>{2, 1, 3}) {
    const std::string var("var" + std::to_string(i));
    eckit::LocalConfiguration conf;
    conf.set("levels", 30);
    vars213SameMeta.push_back(oops::Variable(var, conf));
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
  oops::Variables acb(std::vector<std::string>{"a", "c", "b"});
  oops::Variables b({"b"});
  oops::Variables ba(std::vector<std::string>{"b", "a"});
  oops::Variables de(std::vector<std::string>{"d", "e"});

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
}

// -----------------------------------------------------------------------------

class Variables : public oops::Test {
 public:
  using oops::Test::Test;
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
    ts.emplace_back(CASE("Variables/testVariableConstructorAndEqualsComparison")
      { testVariableConstructorAndEqualsComparison(); });
    ts.emplace_back(CASE("Variables/testPushBack")
      { testPushBack(); });
  }

  void clear() const override {}
};

}  // namespace test

#endif  // TEST_BASE_VARIABLES_H_
