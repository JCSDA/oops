/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_PARAMETERS_H_
#define TEST_UTIL_PARAMETERS_H_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/parameters/ParameterTraitsScalarOrMap.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace test {

enum class Fruit {
  APPLE, ORANGE
};

struct FruitParameterTraitsHelper {
  typedef Fruit EnumType;
  static constexpr char enumTypeName[] = "Fruit";
  static constexpr std::pair<Fruit, const char*> valuesAndNames[] = {
    { Fruit::APPLE, "apple" },
    { Fruit::ORANGE, "orange" }
  };
};

constexpr char FruitParameterTraitsHelper::enumTypeName[];
constexpr std::pair<Fruit, const char*> FruitParameterTraitsHelper::valuesAndNames[];

}  // namespace test

namespace oops {

template <>
struct ParameterTraits<test::Fruit> : public EnumParameterTraits<test::FruitParameterTraitsHelper>
{};

}  // namespace oops

namespace test {

// Classes required by most tests in this file

class RangeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(RangeParameters, Parameters)
 public:
  oops::Parameter<float> minParameter{"min", 0.0f, this};
  oops::Parameter<float> maxParameter{"max", 0.0f, this};
};

class EmbeddedParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(EmbeddedParameters, Parameters)
 public:
  oops::Parameter<int> intParameter{"embedded_int_parameter", 3, this};
  oops::OptionalParameter<util::DateTime> optDateTimeParameter{"opt_embedded_date_time_parameter",
                                                               this};
};

class MyParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(MyParametersBase, Parameters)
 public:
  oops::Parameter<float> floatParameter{"float_parameter", 1.5f, this};
  oops::Parameter<int> intParameter{"int_parameter", 2, this};
  oops::Parameter<bool> boolParameter{"bool_parameter", true, this};
  oops::OptionalParameter<float> optFloatParameter{"opt_float_parameter", this};
  oops::OptionalParameter<util::DateTime> optDateTimeParameter{"opt_date_time_parameter", this};
  oops::OptionalParameter<util::Duration> optDurationParameter{"opt_duration_parameter", this};
  oops::Parameter<Fruit> fruitParameter{"fruit_parameter", Fruit::ORANGE, this};
  oops::Parameter<RangeParameters> rangeParameter{"range_parameter", RangeParameters(), this};
  oops::Parameter<std::vector<int>> intParameters{"int_parameters", {}, this};
  oops::Parameter<std::vector<RangeParameters>> rangeParameters{"range_parameters", {}, this};
  EmbeddedParameters embeddedParameters{this};
};

class MyOptionalParameters : public MyParametersBase {
  OOPS_CONCRETE_PARAMETERS(MyOptionalParameters, MyParametersBase)
 public:
};

class MyOptionalAndRequiredParameters : public MyParametersBase {
  OOPS_CONCRETE_PARAMETERS(MyOptionalAndRequiredParameters, MyParametersBase)
 public:
  oops::RequiredParameter<float> reqFloatParameter{"req_float_parameter", this};
  oops::RequiredParameter<util::Duration> reqDurationParameter{"req_duration_parameter", this};
};

// Class required by tests checking map-valued parameters

class MyMapParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MyMapParameters, Parameters)
 public:
  oops::Parameter<std::map<int, float>> intToFloatMapParameter{
      "int_to_float_map", std::map<int, float>(), this};
  oops::Parameter<std::map<std::string, util::Duration>> stringToDurationMapParameter{
      "string_to_duration_map", std::map<std::string, util::Duration>(), this};

  oops::Parameter<util::ScalarOrMap<int, float>> floatOrIntToFloatMapParameter1{
      "float_or_int_to_float_map_1", util::ScalarOrMap<int, float>(), this};
  oops::Parameter<util::ScalarOrMap<std::string, util::Duration>>
    durationOrStringToDurationMapParameter1{"duration_or_string_to_duration_map_1",
                                            util::ScalarOrMap<std::string, util::Duration>(), this};

  oops::Parameter<util::ScalarOrMap<int, float>> floatOrIntToFloatMapParameter2{
      "float_or_int_to_float_map_2", util::ScalarOrMap<int, float>(), this};
  oops::Parameter<util::ScalarOrMap<std::string, util::Duration>>
    durationOrStringToDurationMapParameter2{"duration_or_string_to_duration_map_2",
                                            util::ScalarOrMap<std::string, util::Duration>(), this};
};

void testDefaultValues() {
  const eckit::LocalConfiguration conf(TestEnvironment::config());

  MyOptionalAndRequiredParameters params;

  EXPECT_EQUAL(params.floatParameter, 1.5f);
  EXPECT_EQUAL(params.floatParameter.value(), 1.5f);
  EXPECT_EQUAL(params.intParameter, 2);
  EXPECT(params.boolParameter);
  EXPECT(params.optFloatParameter.value() == boost::none);
  EXPECT(params.optDateTimeParameter.value() == boost::none);
  EXPECT(params.optDurationParameter.value() == boost::none);
  EXPECT_THROWS_AS(params.reqFloatParameter.value(), boost::bad_optional_access);
  EXPECT_THROWS_AS(params.reqDurationParameter.value(), boost::bad_optional_access);
  EXPECT(params.fruitParameter == Fruit::ORANGE);
  EXPECT(params.rangeParameter.value().minParameter == 0.0f);
  EXPECT(params.rangeParameter.value().maxParameter == 0.0f);
  EXPECT(params.intParameters.value().empty());
  EXPECT(params.rangeParameters.value().empty());
  EXPECT(params.embeddedParameters.intParameter.value() == 3);
  EXPECT(params.embeddedParameters.optDateTimeParameter.value() == boost::none);

  const eckit::LocalConfiguration minimalConf(conf, "minimal");
  params.deserialize(minimalConf);

  EXPECT_EQUAL(params.floatParameter, 1.5f);
  EXPECT_EQUAL(params.floatParameter.value(), 1.5f);
  EXPECT_EQUAL(params.intParameter, 2);
  EXPECT(params.boolParameter);
  EXPECT(params.optFloatParameter.value() == boost::none);
  EXPECT(params.optDateTimeParameter.value() == boost::none);
  EXPECT(params.optDurationParameter.value() == boost::none);
  EXPECT_EQUAL(params.reqFloatParameter, 3.0f);
  EXPECT_EQUAL(params.reqFloatParameter.value(), 3.0f);
  EXPECT_EQUAL(params.reqDurationParameter.value(), util::Duration("PT1H"));
  EXPECT(params.fruitParameter == Fruit::ORANGE);
  EXPECT(params.rangeParameter.value().minParameter == 0.0f);
  EXPECT(params.rangeParameter.value().maxParameter == 0.0f);
  EXPECT(params.intParameters.value().empty());
  EXPECT(params.rangeParameters.value().empty());
  EXPECT(params.embeddedParameters.intParameter.value() == 3);
  EXPECT(params.embeddedParameters.optDateTimeParameter.value() == boost::none);
}

void testCorrectValues() {
  MyOptionalAndRequiredParameters params;
  const eckit::LocalConfiguration fullConf(TestEnvironment::config(), "full");
  params.deserialize(fullConf);

  EXPECT_EQUAL(params.floatParameter, 3.5f);
  EXPECT_EQUAL(params.floatParameter.value(), 3.5f);
  EXPECT_EQUAL(params.intParameter, 4);
  EXPECT(!params.boolParameter);
  EXPECT(params.optFloatParameter.value() != boost::none);
  EXPECT_EQUAL(params.optFloatParameter.value().get(), 5.5f);
  EXPECT(params.optDateTimeParameter.value() != boost::none);
  EXPECT_EQUAL(params.optDateTimeParameter.value().get(), util::DateTime(2010, 2, 3, 4, 5, 6));
  EXPECT(params.optDurationParameter.value() != boost::none);
  EXPECT_EQUAL(params.optDurationParameter.value().get(), util::Duration("PT01H02M03S"));
  EXPECT_EQUAL(params.reqFloatParameter, 6.0f);
  EXPECT_EQUAL(params.reqFloatParameter.value(), 6.0f);
  EXPECT_EQUAL(params.reqDurationParameter.value(), util::Duration("PT06H30M"));
  EXPECT(params.fruitParameter == Fruit::APPLE);
  EXPECT(params.rangeParameter.value().minParameter == 7.0f);
  EXPECT(params.rangeParameter.value().maxParameter == 8.5f);
  EXPECT(params.intParameters.value() == std::vector<int>({1, 2}));
  EXPECT(params.rangeParameters.value().size() == 2);
  EXPECT(params.rangeParameters.value()[0].minParameter == 9.0f);
  EXPECT(params.rangeParameters.value()[0].maxParameter == 10.0f);
  EXPECT(params.rangeParameters.value()[1].minParameter == 11.0f);
  EXPECT(params.rangeParameters.value()[1].maxParameter == 12.0f);
  EXPECT(params.embeddedParameters.intParameter.value() == 13);
  EXPECT(params.embeddedParameters.optDateTimeParameter.value() != boost::none);
  EXPECT_EQUAL(params.embeddedParameters.optDateTimeParameter.value().get(),
               util::DateTime(2010, 3, 4, 5, 6, 7));
}

void testIncorrectValueOfFloatParameter() {
  MyOptionalAndRequiredParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_float_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfOptionalFloatParameter() {
  MyOptionalAndRequiredParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_opt_float_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfOptionalDateTimeParameter() {
  MyOptionalAndRequiredParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_opt_date_time_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfOptionalDurationParameter() {
  MyOptionalAndRequiredParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_opt_duration_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfEnumParameter() {
  MyOptionalAndRequiredParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_fruit_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfIntParameters() {
  MyOptionalAndRequiredParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_int_parameters");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::Exception);
}

void testIncorrectValueOfRangeParameters() {
  MyOptionalAndRequiredParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_range_parameters");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::Exception);
}

void testMissingRequiredFloatParameter() {
  MyOptionalAndRequiredParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "missing_req_float_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testMissingRequiredDurationParameter() {
  MyOptionalAndRequiredParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "missing_req_duration_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

// Tests of parameters storing maps and ScalarOrMaps

void testMapParameters(const MyMapParameters &params) {
  // Map parameters
  EXPECT_EQUAL(params.intToFloatMapParameter.value().at(5),
               1.5f);
  EXPECT_EQUAL(params.intToFloatMapParameter.value().at(7),
               3.0f);
  EXPECT(params.intToFloatMapParameter.value().find(123456) ==
         params.intToFloatMapParameter.value().end());

  EXPECT_EQUAL(params.stringToDurationMapParameter.value().at("day"),
               util::Duration("PT16H"));
  EXPECT_EQUAL(params.stringToDurationMapParameter.value().at("night"),
               util::Duration("PT8H"));
  EXPECT(params.stringToDurationMapParameter.value().find("abcdef") ==
         params.stringToDurationMapParameter.value().end());

  // Scalar-or-map parameters set to maps
  EXPECT_NOT(params.floatOrIntToFloatMapParameter1.value().isScalar());
  EXPECT_EQUAL(params.floatOrIntToFloatMapParameter1.value().at(6),
               2.5f);
  EXPECT_EQUAL(params.floatOrIntToFloatMapParameter1.value().at(8),
               4.0f);
  EXPECT_NOT(params.floatOrIntToFloatMapParameter1.value().contains(1));

  EXPECT_NOT(params.durationOrStringToDurationMapParameter1.value().isScalar());
  EXPECT_EQUAL(params.durationOrStringToDurationMapParameter1.value().at("day"),
               util::Duration("PT14H"));
  EXPECT_EQUAL(params.durationOrStringToDurationMapParameter1.value().at("night"),
               util::Duration("PT10H"));
  EXPECT_NOT(params.durationOrStringToDurationMapParameter1.value().contains("abcdef"));

  // Scalar-or-map parameters set to scalars
  EXPECT(params.floatOrIntToFloatMapParameter2.value().isScalar());
  EXPECT_EQUAL(params.floatOrIntToFloatMapParameter2.value().at(6),
               3.5f);
  EXPECT_EQUAL(params.floatOrIntToFloatMapParameter2.value().at(8),
               3.5f);
  EXPECT(params.floatOrIntToFloatMapParameter2.value().contains(123456));
  EXPECT_EQUAL(params.floatOrIntToFloatMapParameter2.value().at(123456),
               3.5f);

  EXPECT(params.durationOrStringToDurationMapParameter2.value().isScalar());
  EXPECT_EQUAL(params.durationOrStringToDurationMapParameter2.value().at("day"),
               util::Duration("PT12H"));
  EXPECT_EQUAL(params.durationOrStringToDurationMapParameter2.value().at("night"),
               util::Duration("PT12H"));
  EXPECT(params.durationOrStringToDurationMapParameter2.value().contains("abcdef"));
  EXPECT_EQUAL(params.durationOrStringToDurationMapParameter2.value().at("abcdef"),
               util::Duration("PT12H"));
}

void testMapParametersYamlStyleQuotedKeys() {
  MyMapParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "map_parameter_yaml_style_quoted_keys");
  params.deserialize(conf);
  testMapParameters(params);
}

void testMapParametersYamlStyleUnquotedKeys() {
  MyMapParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "map_parameter_yaml_style_unquoted_keys");
  params.deserialize(conf);
  testMapParameters(params);
}

void testMapParametersJsonStyleQuotedKeys() {
  MyMapParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "map_parameter_json_style_quoted_keys");
  params.deserialize(conf);
  testMapParameters(params);
}

void testMapParametersJsonStyleUnquotedKeys() {
  MyMapParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "map_parameter_json_style_unquoted_keys");
  params.deserialize(conf);
  testMapParameters(params);
}

// Tests of special member functions

void expectMatchesFullConf(const MyOptionalAndRequiredParameters &params) {
  EXPECT_EQUAL(params.floatParameter, 3.5f);
  EXPECT_EQUAL(params.rangeParameter.value().minParameter, 7.0f);
  EXPECT_EQUAL(params.rangeParameters.value()[0].minParameter, 9.0f);
  EXPECT(params.embeddedParameters.intParameter.value() == 13);
}

void expectMatchesAlternativeConf(const MyOptionalAndRequiredParameters &params) {
  EXPECT_EQUAL(params.floatParameter, 13.5f);
  EXPECT_EQUAL(params.rangeParameter.value().minParameter, 17.0f);
  EXPECT_EQUAL(params.rangeParameters.value()[0].minParameter, 19.0f);
  EXPECT(params.embeddedParameters.intParameter.value() == 23);
}

void testCopyConstructor() {
  const eckit::LocalConfiguration fullConf(TestEnvironment::config(), "full");
  const eckit::LocalConfiguration alternativeConf(TestEnvironment::config(), "alternative");

  MyOptionalAndRequiredParameters params;
  params.deserialize(fullConf);

  MyOptionalAndRequiredParameters otherParams = params;
  expectMatchesFullConf(otherParams);

  otherParams.deserialize(alternativeConf);
  expectMatchesAlternativeConf(otherParams);

  expectMatchesFullConf(params);

  params.deserialize(alternativeConf);
  expectMatchesAlternativeConf(params);
}

void testMoveConstructor() {
  const eckit::LocalConfiguration fullConf(TestEnvironment::config(), "full");
  const eckit::LocalConfiguration alternativeConf(TestEnvironment::config(), "alternative");

  MyOptionalAndRequiredParameters params;
  params.deserialize(fullConf);

  MyOptionalAndRequiredParameters otherParams = std::move(params);
  expectMatchesFullConf(otherParams);

  otherParams.deserialize(alternativeConf);
  expectMatchesAlternativeConf(otherParams);
}

void testCopyAssignmentOperator() {
  const eckit::LocalConfiguration fullConf(TestEnvironment::config(), "full");
  const eckit::LocalConfiguration alternativeConf(TestEnvironment::config(), "alternative");

  MyOptionalAndRequiredParameters params;
  params.deserialize(fullConf);

  MyOptionalAndRequiredParameters otherParams;
  otherParams = params;
  expectMatchesFullConf(otherParams);

  otherParams.deserialize(alternativeConf);
  expectMatchesAlternativeConf(otherParams);

  expectMatchesFullConf(params);

  params.deserialize(alternativeConf);
  expectMatchesAlternativeConf(params);
}

void testMoveAssignmentOperator() {
  const eckit::LocalConfiguration fullConf(TestEnvironment::config(), "full");
  const eckit::LocalConfiguration alternativeConf(TestEnvironment::config(), "alternative");

  MyOptionalAndRequiredParameters params;
  params.deserialize(fullConf);

  MyOptionalAndRequiredParameters otherParams;
  otherParams = std::move(params);
  expectMatchesFullConf(otherParams);

  otherParams.deserialize(alternativeConf);
  expectMatchesAlternativeConf(otherParams);
}

class Parameters : public oops::Test {
 private:
  std::string testid() const override {return "test::Parameters";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("util/Parameters/defaultValues") {
                      testDefaultValues();
                    });
    ts.emplace_back(CASE("util/Parameters/correctValues") {
                      testCorrectValues();
                    });
    ts.emplace_back(CASE("util/Parameters/incorrectValueOfFloatParameter") {
                      testIncorrectValueOfFloatParameter();
                    });
    ts.emplace_back(CASE("util/Parameters/incorrectValueOfOptionalFloatParameter") {
                      testIncorrectValueOfOptionalFloatParameter();
                    });
    // Conversion from string to DateTime or Duration calls abort() on failure,
    // so we can't test these cases. Leaving them commented-out in case this changes in future.
    //
    // ts.emplace_back(CASE("util/Parameters/incorrectValueOfOptionalDateTimeParameter") {
    //                   testIncorrectValueOfOptionalDateTimeParameter();
    //                 });
    // ts.emplace_back(CASE("util/Parameters/incorrectValueOfOptionalDurationParameter") {
    //                   testIncorrectValueOfOptionalDurationParameter();
    //                 });
    ts.emplace_back(CASE("util/Parameters/incorrectValueOfEnumParameter") {
                      testIncorrectValueOfEnumParameter();
                    });
    ts.emplace_back(CASE("util/Parameters/testIncorrectValueOfIntParameters") {
                      testIncorrectValueOfIntParameters();
                    });
    ts.emplace_back(CASE("util/Parameters/testIncorrectValueOfRangeParameters") {
                      testIncorrectValueOfRangeParameters();
                    });
    ts.emplace_back(CASE("util/Parameters/testMissingRequiredFloatParameter") {
                      testMissingRequiredFloatParameter();
                    });
    ts.emplace_back(CASE("util/Parameters/testMissingRequiredDurationParameter") {
                      testMissingRequiredDurationParameter();
                    });

    ts.emplace_back(CASE("util/Parameters/testCopyConstructor") {
                      testCopyConstructor();
                    });
    ts.emplace_back(CASE("util/Parameters/testMoveConstructor") {
                      testMoveConstructor();
                    });
    ts.emplace_back(CASE("util/Parameters/testCopyAssignmentOperator") {
                      testCopyAssignmentOperator();
                    });
    ts.emplace_back(CASE("util/Parameters/testMoveAssignmentOperator") {
                      testMoveAssignmentOperator();
                    });

    // Test fails because of a bug in the eckit YAML parser
    // ts.emplace_back(CASE("util/Parameters/mapParametersYamlStyleQuotedKeys") {
    //                   testMapParametersYamlStyleQuotedKeys();
    //                 });
    // Test fails because of a bug in the eckit YAML parser
    // ts.emplace_back(CASE("util/Parameters/mapParametersYamlStyleUnquotedKeys") {
    //                   testMapParametersYamlStyleUnquotedKeys();
    //                 });
    ts.emplace_back(CASE("util/Parameters/mapParametersJsonStyleQuotedKeys") {
                      testMapParametersJsonStyleQuotedKeys();
                    });
    // Test fails because of a bug in the eckit YAML parser
    // ts.emplace_back(CASE("util/Parameters/mapParametersJsonStyleUnquotedKeys") {
    //                   testMapParametersJsonStyleUnquotedKeys();
    //                 });
  }
};

}  // namespace test

#endif  // TEST_UTIL_PARAMETERS_H_
