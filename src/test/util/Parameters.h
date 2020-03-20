/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_PARAMETERS_H_
#define TEST_UTIL_PARAMETERS_H_

#include <map>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/parameters/ParameterTraitsScalarOrMap.h"

namespace test {

enum class Fruit {
  APPLE, ORANGE
};

}  // namespace test

namespace oops {

template <>
struct ParameterTraits<test::Fruit> {
  static boost::optional<test::Fruit> get(const eckit::Configuration &config,
                                                const std::string& name) {
    std::string value;
    if (config.get(name, value)) {
      if (value == "apple")
        return test::Fruit::APPLE;
      if (value == "orange")
        return test::Fruit::ORANGE;
      throw eckit::BadParameter("Bad conversion from std::string '" + value + "' to Fruit",
                                Here());
    } else {
      return boost::none;
    }
  }
};

}  // namespace oops

namespace test {

class RangeParameters : public oops::Parameters {
 public:
  oops::Parameter<float> minParameter{"min", 0.0f, this};
  oops::Parameter<float> maxParameter{"max", 0.0f, this};
};

class MyParameters : public oops::Parameters {
 public:
  oops::Parameter<float> floatParameter{"float_parameter", 1.5f, this};
  oops::Parameter<int> intParameter{"int_parameter", 2, this};
  oops::Parameter<bool> boolParameter{"bool_parameter", true, this};
  oops::OptionalParameter<float> optFloatParameter{"opt_float_parameter", this};
  oops::OptionalParameter<util::DateTime> optDateTimeParameter{"opt_date_time_parameter", this};
  oops::OptionalParameter<util::Duration> optDurationParameter{"opt_duration_parameter", this};
  oops::Parameter<Fruit> fruitParameter{"fruit_parameter", Fruit::ORANGE, this};
  oops::Parameter<RangeParameters> rangeParameter{"range_parameter", {}, this};
  oops::Parameter<std::vector<int>> intParameters{"int_parameters", {}, this};
  oops::Parameter<std::vector<RangeParameters>> rangeParameters{"range_parameters", {}, this};
};

class MyMapParameters : public oops::Parameters {
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

  MyParameters params;

  EXPECT_EQUAL(params.floatParameter, 1.5f);
  EXPECT_EQUAL(params.intParameter, 2);
  EXPECT(params.boolParameter);
  EXPECT(params.optFloatParameter.value() == boost::none);
  EXPECT(params.optDateTimeParameter.value() == boost::none);
  EXPECT(params.optDurationParameter.value() == boost::none);
  EXPECT(params.fruitParameter == Fruit::ORANGE);
  EXPECT(params.rangeParameter.value().minParameter == 0.0f);
  EXPECT(params.rangeParameter.value().maxParameter == 0.0f);
  EXPECT(params.intParameters.value().empty());
  EXPECT(params.rangeParameters.value().empty());

  const eckit::LocalConfiguration emptyConf(conf, "empty");
  params.deserialize(emptyConf);

  EXPECT_EQUAL(params.floatParameter, 1.5f);
  EXPECT_EQUAL(params.intParameter, 2);
  EXPECT(params.boolParameter);
  EXPECT(params.optFloatParameter.value() == boost::none);
  EXPECT(params.optDateTimeParameter.value() == boost::none);
  EXPECT(params.optDurationParameter.value() == boost::none);
  EXPECT(params.fruitParameter == Fruit::ORANGE);
  EXPECT(params.rangeParameter.value().minParameter == 0.0f);
  EXPECT(params.rangeParameter.value().maxParameter == 0.0f);
  EXPECT(params.intParameters.value().empty());
  EXPECT(params.rangeParameters.value().empty());
}

void testCorrectValues() {
  MyParameters params;
  const eckit::LocalConfiguration fullConf(TestEnvironment::config(), "full");
  params.deserialize(fullConf);

  EXPECT_EQUAL(params.floatParameter, 3.5f);
  EXPECT_EQUAL(params.intParameter, 4);
  EXPECT(!params.boolParameter);
  EXPECT(params.optFloatParameter.value() != boost::none);
  EXPECT_EQUAL(params.optFloatParameter.value().get(), 5.5f);
  EXPECT(params.optDateTimeParameter.value() != boost::none);
  EXPECT_EQUAL(params.optDateTimeParameter.value().get(), util::DateTime(2010, 2, 3, 4, 5, 6));
  EXPECT(params.optDurationParameter.value() != boost::none);
  EXPECT_EQUAL(params.optDurationParameter.value().get(), util::Duration("PT01H02M03S"));
  EXPECT(params.fruitParameter == Fruit::APPLE);
  EXPECT(params.rangeParameter.value().minParameter == 7.0f);
  EXPECT(params.rangeParameter.value().maxParameter == 8.5f);
  EXPECT(params.intParameters.value() == std::vector<int>({1, 2}));
  EXPECT(params.rangeParameters.value().size() == 2);
  EXPECT(params.rangeParameters.value()[0].minParameter == 9.0f);
  EXPECT(params.rangeParameters.value()[0].maxParameter == 10.0f);
  EXPECT(params.rangeParameters.value()[1].minParameter == 11.0f);
  EXPECT(params.rangeParameters.value()[1].maxParameter == 12.0f);
}

void testIncorrectValueOfFloatParameter() {
  MyParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_float_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfOptionalFloatParameter() {
  MyParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_opt_float_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfOptionalDateTimeParameter() {
  MyParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_opt_date_time_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfOptionalDurationParameter() {
  MyParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_opt_duration_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfEnumParameter() {
  MyParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_fruit_parameter");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::BadParameter);
}

void testIncorrectValueOfIntParameters() {
  MyParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_int_parameters");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::Exception);
}

void testIncorrectValueOfRangeParameters() {
  MyParameters params;
  const eckit::LocalConfiguration conf(TestEnvironment::config(),
                                       "error_in_range_parameters");
  EXPECT_THROWS_AS(params.deserialize(conf), eckit::Exception);
}

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
