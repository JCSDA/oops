/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_COMPOSITEPATH_H_
#define TEST_UTIL_COMPOSITEPATH_H_

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/CompositePath.h"
#include "oops/util/Expect.h"

namespace test {

CASE("util/CompositePath/compositePath") {
  util::CompositePath path('\\');
  EXPECT_EQUAL(path.path(), "\\");
  {
    util::PathComponent abc(path, "abc");
    EXPECT_EQUAL(path.path(), "\\abc");
    {
      util::PathComponent abc(path, "def");
      EXPECT_EQUAL(path.path(), "\\abc\\def");
    }
    EXPECT_EQUAL(path.path(), "\\abc");
    {
      util::PathComponent abc(path, "ghij");
      EXPECT_EQUAL(path.path(), "\\abc\\ghij");
    }
    EXPECT_EQUAL(path.path(), "\\abc");
  }
  EXPECT_EQUAL(path.path(), "\\");
  {
    util::PathComponent abc(path, "klm");
    EXPECT_EQUAL(path.path(), "\\klm");
  }
  EXPECT_EQUAL(path.path(), "\\");
}

class CompositePath : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::CompositePath";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_COMPOSITEPATH_H_
