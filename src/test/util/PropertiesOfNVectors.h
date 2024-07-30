/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_PROPERTIESOFNVECTORS_H_
#define TEST_UTIL_PROPERTIESOFNVECTORS_H_

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/PropertiesOfNVectors.h"

namespace test {

  void testSizes()
  {
    // Vectors used in testing
    const std::vector <size_t> vEmpty {};
    const std::vector <int> v1 {1};
    const std::vector <float> v2 {1.0, 2.0};
    const std::vector <std::string> v3 {"1", "2", "3"};

    // Tests with one vector
    EXPECT_EQUAL(oops::listOfVectorSizes(vEmpty), "0");
    EXPECT_EQUAL(oops::listOfVectorSizes(v1), "1");
    EXPECT_EQUAL(oops::listOfVectorSizes(v2), "2");
    EXPECT_EQUAL(oops::listOfVectorSizes(v3), "3");

    // Tests with more than one vector
    EXPECT_EQUAL(oops::listOfVectorSizes(vEmpty, v1), "0, 1");
    EXPECT_EQUAL(oops::listOfVectorSizes(vEmpty, v1, v2), "0, 1, 2");
    EXPECT_EQUAL(oops::listOfVectorSizes(vEmpty, v1, v2, v3), "0, 1, 2, 3");
  }

  CASE("util/PropertiesOfNVectors/sizes") {
    testSizes();
  }

class PropertiesOfNVectors : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::PropertiesOfNVectors";}
  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_PROPERTIESOFNVECTORS_H_
