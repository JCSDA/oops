/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_SCALARORMAP_H_
#define TEST_UTIL_SCALARORMAP_H_

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/ScalarOrMap.h"

namespace test {

CASE("util/ScalarOrMap/storingMap") {
  const eckit::LocalConfiguration conf(TestEnvironment::config());
  typedef util::ScalarOrMap<int, float> MapType;
  MapType map1(std::map<int, float>{{1, 3.14f}, {2, 1.5f}});

  EXPECT_NOT(const_cast<const MapType&>(map1).isScalar());
  EXPECT(const_cast<const MapType&>(map1).contains(1));
  EXPECT_NOT(const_cast<const MapType&>(map1).contains(3));
  EXPECT_EQUAL(const_cast<const MapType&>(map1).at(1), 3.14f);
  EXPECT_EQUAL(const_cast<const MapType&>(map1).at(2), 1.5f);
  EXPECT_THROWS(const_cast<const MapType&>(map1).at(3));

  std::vector<int> keys;
  std::vector<float> values;
  for (const std::pair<const int, float> &kv : map1) {
    keys.push_back(kv.first);
    values.push_back(kv.second);
  }
  EXPECT_EQUAL(keys, (std::vector<int>{1, 2}));
  EXPECT_EQUAL(values, (std::vector<float>{3.14f, 1.5f}));

  keys.clear();
  values.clear();
  for (const std::pair<const int, float> &kv : const_cast<const MapType&>(map1)) {
    keys.push_back(kv.first);
    values.push_back(kv.second);
  }
  EXPECT_EQUAL(keys, (std::vector<int>{1, 2}));
  EXPECT_EQUAL(values, (std::vector<float>{3.14f, 1.5f}));

  map1 = 5.5f;

  EXPECT(const_cast<const MapType&>(map1).isScalar());
  EXPECT_EQUAL(map1.at(1), 5.5f);
}

CASE("util/ScalarOrMap/storingScalar") {
  const eckit::LocalConfiguration conf(TestEnvironment::config());
  typedef util::ScalarOrMap<int, float> MapType;
  MapType map1(2.5f);

  EXPECT(const_cast<const MapType&>(map1).isScalar());
  EXPECT(const_cast<const MapType&>(map1).contains(1));
  EXPECT(const_cast<const MapType&>(map1).contains(3));
  EXPECT_EQUAL(const_cast<const MapType&>(map1).at(1), 2.5f);
  EXPECT_EQUAL(const_cast<const MapType&>(map1).at(2), 2.5f);

  EXPECT_THROWS(map1.begin());
  EXPECT_THROWS(map1.end());
  EXPECT_THROWS(const_cast<const MapType&>(map1).begin());
  EXPECT_THROWS(const_cast<const MapType&>(map1).end());

  map1 = std::map<int, float>{{1, 3.14f}, {2, 1.5f}};

  EXPECT_NOT(const_cast<const MapType&>(map1).isScalar());
  EXPECT_EQUAL(map1.at(1), 3.14f);
}

class ScalarOrMap : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::ScalarOrMap";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_SCALARORMAP_H_
