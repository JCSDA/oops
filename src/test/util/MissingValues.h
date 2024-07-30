/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_MISSINGVALUES_H_
#define TEST_UTIL_MISSINGVALUES_H_

#include <string>

#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/missingValues.h"

namespace test {

template <typename T>
void testMissingValues()
{
  const T missing = util::missingValue<T>();
  (void)missing;  // silence an unused-variable warning
}

CASE("util/MissingValues/float") {
  testMissingValues<float>();
}

CASE("util/MissingValues/double") {
  testMissingValues<double>();
}

CASE("util/MissingValues/int16_t") {
  testMissingValues<int16_t>();
}

CASE("util/MissingValues/int32_t") {
  testMissingValues<int32_t>();
}

CASE("util/MissingValues/int64_t") {
  testMissingValues<int64_t>();
}

CASE("util/MissingValues/bool") {
  testMissingValues<bool>();
}

CASE("util/MissingValues/char") {
  testMissingValues<char>();
}

CASE("util/MissingValues/std::string") {
  testMissingValues<std::string>();
}

CASE("util/MissingValues/DateTime") {
  testMissingValues<util::DateTime>();
}

class MissingValues : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::MissingValues";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_MISSINGVALUES_H_
