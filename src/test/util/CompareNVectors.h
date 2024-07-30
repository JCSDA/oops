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
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"

namespace test {

  void testEmpty()
  {
    // Vectors used in testing
    const std::vector <size_t> vEmpty {};
    const std::vector <int> v1 {1, 2, 3};
    const std::vector <float> v2 {1.0, 2.0, 3.0};
    const std::vector <std::string> v3 {"1", "2", "3"};

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
    const std::vector <int> v1a {1};
    const std::vector <float> v1b {1.0};
    const std::vector <std::string> v1c {"1"};
    const std::vector <int> v2a {1, 2};
    const std::vector <float> v2b {1.0, 2.0};
    const std::vector <std::string> v2c {"1", "2"};

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

  void testExpectedSize()
  {
    // Vectors used in testing
    const std::vector <size_t> vEmpty {};
    const std::vector <int> v1a {1};
    const std::vector <float> v1b {1.0};
    const std::vector <std::string> v1c {"1"};
    const std::vector <int> v2a {1, 2};
    const std::vector <float> v2b {1.0, 2.0};
    const std::vector <std::string> v2c {"1", "2"};

    // Expected vector size (N)
    const size_t N = 2;

    // Vectors of length 2
    EXPECT(oops::allVectorsExpectedSize(N, v2a));
    EXPECT(oops::allVectorsExpectedSize(N, v2a, v2b));
    EXPECT(oops::allVectorsExpectedSize(N, v2a, v2b, v2c));

    // Empty vector
    EXPECT_NOT(oops::allVectorsExpectedSize(N, vEmpty));

    // Vectors of length 1
    EXPECT_NOT(oops::allVectorsExpectedSize(N, v1a));
    EXPECT_NOT(oops::allVectorsExpectedSize(N, v1a, v1b));
    EXPECT_NOT(oops::allVectorsExpectedSize(N, v1a, v1c));

    // Vectors of different lengths
    EXPECT_NOT(oops::allVectorsExpectedSize(N, v2a, vEmpty));
    EXPECT_NOT(oops::allVectorsExpectedSize(N, v2a, v1a, vEmpty));
    EXPECT_NOT(oops::allVectorsExpectedSize(N, v2a, vEmpty, vEmpty));
    EXPECT_NOT(oops::allVectorsExpectedSize(N, vEmpty, v2a, vEmpty, v2b, vEmpty, v2c));
  }

  void testNonEmptyExpectedSize()
  {
    // Vectors used in testing
    const std::vector <size_t> vEmpty {};
    const std::vector <int> v1a {1};
    const std::vector <float> v1b {1.0};
    const std::vector <std::string> v1c {"1"};
    const std::vector <int> v2a {1, 2};
    const std::vector <float> v2b {1.0, 2.0};
    const std::vector <std::string> v2c {"1", "2"};

    // Expected vector size (N)
    const size_t N = 2;

    // Vectors of length 2
    EXPECT(oops::allNonEmptyVectorsExpectedSize(N, v2a));
    EXPECT(oops::allNonEmptyVectorsExpectedSize(N, v2a, v2b));
    EXPECT(oops::allNonEmptyVectorsExpectedSize(N, v2a, v2b, v2c));

    // Empty vector
    EXPECT(oops::allNonEmptyVectorsExpectedSize(N, vEmpty));

    // Vectors either of length 2 or empty
    EXPECT(oops::allNonEmptyVectorsExpectedSize(N, v2a, vEmpty));
    EXPECT(oops::allNonEmptyVectorsExpectedSize(N, v2a, vEmpty, vEmpty));
    EXPECT(oops::allNonEmptyVectorsExpectedSize(N, vEmpty, v2a, vEmpty, v2b, vEmpty, v2c));

    // Vectors of length 1
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v1a));
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v1a, v1b));
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v1a, v1c));

    // Vectors either of length 1 or empty
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v1a, vEmpty));
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v1a, vEmpty, vEmpty));
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v1a, vEmpty, v1c, vEmpty));

    // Vectors of different lengths
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v1a, v2a));
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v2a, v1a));
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v1a, v2a, vEmpty));
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v2a, v1a, vEmpty));
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v1a, v2a, v2b));
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v2a, v1a, v2b));
    EXPECT_NOT(oops::allNonEmptyVectorsExpectedSize(N, v2a, vEmpty, v2b, vEmpty, v1a));
  }

  void testSameNonZeroSize()
  {
    // Vectors used in testing
    const std::vector <size_t> vEmpty1 {};
    const std::vector <std::string> vEmpty2 {};
    const std::vector <int> v1a {1};
    const std::vector <float> v1b {1.0};
    const std::vector <std::string> v1c {"1"};
    const std::vector <int> v2a {1, 2};
    const std::vector <float> v2b {1.0, 2.0};
    const std::vector <std::string> v2c {"1", "2"};

    // Empty vectors
    EXPECT_NOT(oops::allVectorsSameNonZeroSize(vEmpty1));
    EXPECT_NOT(oops::allVectorsSameNonZeroSize(vEmpty1, vEmpty2));

    // Vectors of length 1
    EXPECT(oops::allVectorsSameNonZeroSize(v1a));
    EXPECT(oops::allVectorsSameNonZeroSize(v1a, v1b));
    EXPECT(oops::allVectorsSameNonZeroSize(v1a, v1b, v1c));

    // Vectors of length 2
    EXPECT(oops::allVectorsSameNonZeroSize(v2a));
    EXPECT(oops::allVectorsSameNonZeroSize(v2a, v2b));
    EXPECT(oops::allVectorsSameNonZeroSize(v2a, v2b, v2c));

    // Two vectors of different lengths (one empty)
    EXPECT_NOT(oops::allVectorsSameNonZeroSize(vEmpty1, v1a));

    // Two non-empty vectors of different lengths
    EXPECT_NOT(oops::allVectorsSameNonZeroSize(v1a, v2a));

    // Three vectors of different lengths (one empty)
    EXPECT_NOT(oops::allVectorsSameNonZeroSize(vEmpty1, v1a, v1b));
    EXPECT_NOT(oops::allVectorsSameNonZeroSize(v1a, vEmpty1, v1b));
    EXPECT_NOT(oops::allVectorsSameNonZeroSize(v1a, v1b, vEmpty1));

    // Four non-empty vectors of different lengths, various orders
    EXPECT_NOT(oops::allVectorsSameNonZeroSize(v1a, v2a, v2b, v2c));
    EXPECT_NOT(oops::allVectorsSameNonZeroSize(v2a, v1a, v2b, v2c));
    EXPECT_NOT(oops::allVectorsSameNonZeroSize(v2a, v2b, v1a, v2c));
    EXPECT_NOT(oops::allVectorsSameNonZeroSize(v2a, v2b, v2c, v1a));
  }

  CASE("util/CompareNVectors/empty") {
    testEmpty();
  }

  CASE("util/CompareNVectors/sameSize") {
    testSameSize();
  }

  CASE("util/CompareNVectors/expectedSize") {
    testExpectedSize();
  }

  CASE("util/CompareNVectors/nonEmptyExpectedSize") {
    testNonEmptyExpectedSize();
  }

  CASE("util/CompareNVectors/sameNonZeroSize") {
    testSameNonZeroSize();
  }

class CompareNVectors : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::CompareNVectors";}
  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_COMPARENVECTORS_H_
