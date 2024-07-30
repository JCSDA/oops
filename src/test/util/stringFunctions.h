/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_STRINGFUNCTIONS_H_
#define TEST_UTIL_STRINGFUNCTIONS_H_

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/stringFunctions.h"

namespace test {

CASE("util/stringFunctions/join/emptyVector") {
  std::vector<int> v;
  const std::string result = util::stringfunctions::join(
        ", ", v.begin(), v.end(), [](int n){ return std::to_string(n); });
  EXPECT_EQUAL(result, "");
}

CASE("util/stringFunctions/join/singleElementVector") {
  std::vector<int> v{1};
  const std::string result = util::stringfunctions::join(
        ", ", v.begin(), v.end(), [](int n){ return std::to_string(n); });
  EXPECT_EQUAL(result, "1");
}

CASE("util/stringFunctions/join/twoElementVector") {
  std::vector<int> v{1, 2};
  const std::string result = util::stringfunctions::join(
        ", ", v.begin(), v.end(), [](int n){ return std::to_string(n); });
  EXPECT_EQUAL(result, "1, 2");
}

class StringFunctions : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::stringFunctions";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_STRINGFUNCTIONS_H_
