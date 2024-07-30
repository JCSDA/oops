/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_RANGE_H_
#define TEST_UTIL_RANGE_H_

#include <string>

#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/Range.h"

namespace test {

CASE("util/Range/equality") {
  util::Range<int> a{2, 9};
  util::Range<int> b{0, 9};
  util::Range<int> c{2, 5};
  util::Range<int> d{2, 9};

  EXPECT_NOT(a == b);
  EXPECT_NOT(a == c);
  EXPECT(a == d);
}

CASE("util/Range/output") {
  std::stringstream str;
  str << util::Range<int>{2, 9};

  std::string output = str.str();
  EXPECT_EQUAL(output, "[2, 9)");
}

class Range : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::Range";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_RANGE_H_
