/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_ALGORITHMS_H_
#define TEST_UTIL_ALGORITHMS_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/algorithms.h"
#include "oops/util/Expect.h"

namespace test {

CASE("util/Algorithms/transformVector") {
  const std::vector<int> input1{1, 2, 3};
  const std::vector<std::string> output1 =
      util::transformVector(input1, [](int x) { return std::to_string(x); });
  const std::vector<std::string> expectedOutput1{"1", "2", "3"};
  EXPECT_EQUAL(output1, expectedOutput1);

  const std::vector<int> input2;
  const std::vector<std::string> output2 =
      util::transformVector(input2, [](int x) { return std::to_string(x); });
  EXPECT(output2.empty());
}

class Algorithms : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override { return "test::algorithms"; }
  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_ALGORITHMS_H_
