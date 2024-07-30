/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_TYPETRAITS_H_
#define TEST_UTIL_TYPETRAITS_H_

#include <string>
#include <vector>

#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/TypeTraits.h"

namespace test {

CASE("util/TypeTraits/any_is_same") {
  // The extra parentheses are necessary because EXPECT is a macro and its argument contains commas.
  EXPECT((util::any_is_same<int, int>::value));
  EXPECT((util::any_is_same<std::vector<float>, std::vector<float>>::value));
  EXPECT((util::any_is_same<int, int, std::vector<float>>::value));
  EXPECT((util::any_is_same<std::vector<float>, int, std::vector<float>>::value));
  EXPECT((util::any_is_same<std::vector<float>, int, std::vector<int>, std::vector<float>>::value));
  EXPECT_NOT((util::any_is_same<int, const int>::value));
  EXPECT_NOT((util::any_is_same<int, float>::value));
  EXPECT_NOT((util::any_is_same<std::vector<float>, int, std::vector<int>>::value));
}

class TypeTraits : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::TypeTraits";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_TYPETRAITS_H_
