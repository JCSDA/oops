/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_ASSOCIATIVECONTAINERS_H_
#define TEST_UTIL_ASSOCIATIVECONTAINERS_H_

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Expect.h"

namespace test {

CASE("util/AssociativeContainers/intMapContains") {
  std::map<int, float> map{{1, 2.0f}, {3, 4.0f}};

  EXPECT(oops::contains(map, 1));
  EXPECT(oops::contains(map, 3));
  EXPECT_NOT(oops::contains(map, 0));
  EXPECT_NOT(oops::contains(map, 2));
}

CASE("util/AssociativeContainers/stringMapContains") {
  std::map<std::string, float> map{{"abc", 2.0f}, {"def", 4.0f}};

  // This tests if the "key" argument can be of a type convertible to the map key, but not identical
  // to it (const char* vs. std::string).
  EXPECT(oops::contains(map, "abc"));
  EXPECT(oops::contains(map, "def"));
  EXPECT_NOT(oops::contains(map, ""));
  EXPECT_NOT(oops::contains(map, "xxx"));
}

CASE("util/AssociativeContainers/intSetContains") {
  std::set<int> set{1, 3};

  EXPECT(oops::contains(set, 1));
  EXPECT(oops::contains(set, 3));
  EXPECT_NOT(oops::contains(set, 0));
  EXPECT_NOT(oops::contains(set, 2));
}

CASE("util/AssociativeContainers/stringSetContains") {
  std::set<std::string> set{"abc", "def"};

  EXPECT(oops::contains(set, "abc"));
  EXPECT(oops::contains(set, "def"));
  EXPECT_NOT(oops::contains(set, ""));
  EXPECT_NOT(oops::contains(set, "xxx"));
}

CASE("util/AssociativeContainers/keysOfEmptyMap") {
  std::map<std::string, float> map;
  std::vector<std::string> expectedKeys;
  EXPECT_EQUAL(oops::keys(map), expectedKeys);
}

CASE("util/AssociativeContainers/keysOfNonemptyMap") {
  std::map<std::string, float> map{{"abc", 2.0f}, {"def", 4.0f}};
  std::vector<std::string> expectedKeys{"abc", "def"};
  EXPECT_EQUAL(oops::keys(map), expectedKeys);
}

class AssociativeContainers : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::AssociativeContainers";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_ASSOCIATIVECONTAINERS_H_
