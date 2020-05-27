/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_COMPARENVECTORS_H_
#define TEST_UTIL_COMPARENVECTORS_H_

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"

namespace test {

  void testEmpty()
  {
    // Vectors used in testing
    std::vector <size_t> vEmpty {};
    std::vector <int> v1 {1, 2, 3};
    std::vector <float> v2 {1.0, 2.0, 3.0};
    std::vector <std::string> v3 {"1", "2", "3"};

    // Tests with one empty vector
    EXPECT(oops::anyVectorEmpty(vEmpty));
    EXPECT(oops::anyVectorEmpty(vEmpty, v1));
    EXPECT(oops::anyVectorEmpty(vEmpty, v1, v2));
    EXPECT(oops::anyVectorEmpty(vEmpty, v1, v2, v3));

    // Tests with one empty vector and different ordering
    EXPECT(oops::anyVectorEmpty(v1, vEmpty));
    EXPECT(oops::anyVectorEmpty(v1, vEmpty, v2));
    EXPECT(oops::anyVectorEmpty(v1, v2, vEmpty));
    EXPECT(oops::anyVectorEmpty(v1, vEmpty, v2, v3));
    EXPECT(oops::anyVectorEmpty(v1, v2, vEmpty, v3));
    EXPECT(oops::anyVectorEmpty(v1, v2, v3, vEmpty));

    // Tests with two empty vectors in various orders
    EXPECT(oops::anyVectorEmpty(vEmpty, vEmpty));
    EXPECT(oops::anyVectorEmpty(vEmpty, vEmpty, v1));
    EXPECT(oops::anyVectorEmpty(vEmpty, v1, vEmpty));
    EXPECT(oops::anyVectorEmpty(v1, vEmpty, vEmpty));

    // Tests with no empty vectors
    EXPECT_NOT(oops::anyVectorEmpty(v1));
    EXPECT_NOT(oops::anyVectorEmpty(v1, v2));
    EXPECT_NOT(oops::anyVectorEmpty(v1, v2, v3));
  }

  void testSameSize()
  {
    // Vectors used in testing
    std::vector <int> v1a {1};
    std::vector <float> v1b {1.0};
    std::vector <std::string> v1c {"1"};
    std::vector <int> v2a {1, 2};
    std::vector <float> v2b {1.0, 2.0};
    std::vector <std::string> v2c {"1", "2"};

    // Vectors of length 1
    EXPECT(oops::allVectorsSameSize(v1a));
    EXPECT(oops::allVectorsSameSize(v1a, v1b));
    EXPECT(oops::allVectorsSameSize(v1a, v1b, v1c));

    // Vectors of length 2
    EXPECT(oops::allVectorsSameSize(v2a));
    EXPECT(oops::allVectorsSameSize(v2a, v2b));
    EXPECT(oops::allVectorsSameSize(v2a, v2b, v2c));

    // Two vectors of different lengths
    EXPECT_NOT(oops::allVectorsSameSize(v1a, v2a));

    // Four vectors of different lengths, various orders
    EXPECT_NOT(oops::allVectorsSameSize(v1a, v2a, v2b, v2c));
    EXPECT_NOT(oops::allVectorsSameSize(v2a, v1a, v2b, v2c));
    EXPECT_NOT(oops::allVectorsSameSize(v2a, v2b, v1a, v2c));
    EXPECT_NOT(oops::allVectorsSameSize(v2a, v2b, v2c, v1a));
  }

  CASE("util/CompareNVectors/empty") {
    testEmpty();
  }

  CASE("util/CompareNVectors/sameSize") {
    testSameSize();
  }

class CompareNVectors : public oops::Test {
 private:
  std::string testid() const override {return "test::CompareNVectors";}
  void register_tests() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_COMPARENVECTORS_H_
