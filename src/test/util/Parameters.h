/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_PARAMETERS_H_
#define TEST_UTIL_PARAMETERS_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include <boost/make_shared.hpp>

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

class MyParameters : public oops::Parameters {
 public:
  oops::Parameter<float> floatParameter{"float_parameter", 1.5f, this};
  oops::Parameter<int> intParameter{"int_parameter", 2, this};
  oops::Parameter<bool> boolParameter{"bool_parameter", true, this};
  oops::OptionalParameter<float> optFloatParameter{"opt_float_parameter", this};
  oops::OptionalParameter<util::DateTime> optDateTimeParameter{"opt_date_time_parameter", this};
  oops::OptionalParameter<util::Duration> optDurationParameter{"opt_duration_parameter", this};
  oops::Parameter<Fruit> fruitParameter{"fruit_parameter", Fruit::ORANGE, this};
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

  const eckit::LocalConfiguration emptyConf(conf, "empty");
  params.deserialize(emptyConf);

  EXPECT_EQUAL(params.floatParameter, 1.5f);
  EXPECT_EQUAL(params.intParameter, 2);
  EXPECT(params.boolParameter);
  EXPECT(params.optFloatParameter.value() == boost::none);
  EXPECT(params.optDateTimeParameter.value() == boost::none);
  EXPECT(params.optDurationParameter.value() == boost::none);
  EXPECT(params.fruitParameter == Fruit::ORANGE);
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
  }
};

}  // namespace test

#endif  // TEST_UTIL_PARAMETERS_H_
