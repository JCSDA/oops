/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_LOCALENVIRONMENT_H_
#define TEST_UTIL_LOCALENVIRONMENT_H_

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/LocalEnvironment.h"

namespace test {

CASE("util/LocalEnvironment/set") {
  setenv("SOME_VARIABLE", "abcdef", 1 /*replace?*/);
  setenv("ANOTHER_VARIABLE", "XYZ", 1 /*replace?*/);

  {
    util::LocalEnvironment localEnv;

    localEnv.set("SOME_VARIABLE", "1234");
    EXPECT_EQUAL(std::string(::getenv("SOME_VARIABLE")), "1234");
    localEnv.set("NEW_VARIABLE", "PQR");
    EXPECT_EQUAL(std::string(::getenv("NEW_VARIABLE")), "PQR");

    localEnv.set("SOME_VARIABLE", "5678");
    EXPECT_EQUAL(std::string(::getenv("SOME_VARIABLE")), "5678");
    localEnv.set("NEW_VARIABLE", "ijk");
    EXPECT_EQUAL(std::string(::getenv("NEW_VARIABLE")), "ijk");

    // Sanity check: ANOTHER_VARIABLE should still be set to its original value
    EXPECT_EQUAL(std::string(::getenv("ANOTHER_VARIABLE")), "XYZ");
  }

  // SOME_VARIABLE should now be restored to its original value and NEW_VARIABLE unset
  EXPECT_EQUAL(std::string(::getenv("SOME_VARIABLE")), "abcdef");
  EXPECT(::getenv("NEW_VARIABLE") == nullptr);

  // Sanity check: ANOTHER_VARIABLE should still be set to its original value
  EXPECT_EQUAL(std::string(::getenv("ANOTHER_VARIABLE")), "XYZ");
}

CASE("util/LocalEnvironment/exceptionSafety") {
  setenv("SOME_VARIABLE", "abcdef", 1 /*replace?*/);
  setenv("ANOTHER_VARIABLE", "XYZ", 1 /*replace?*/);

  try {
    util::LocalEnvironment localEnv;

    localEnv.set("SOME_VARIABLE", "1234");
    EXPECT_EQUAL(std::string(::getenv("SOME_VARIABLE")), "1234");
    localEnv.set("NEW_VARIABLE", "PQR");
    EXPECT_EQUAL(std::string(::getenv("NEW_VARIABLE")), "PQR");

    throw std::runtime_error("test");
  } catch (std::runtime_error &) {
    // SOME_VARIABLE should now be restored to its original value and NEW_VARIABLE unset
    EXPECT_EQUAL(std::string(::getenv("SOME_VARIABLE")), "abcdef");
    EXPECT(::getenv("NEW_VARIABLE") == nullptr);

    // Sanity check: ANOTHER_VARIABLE should still be set to its original value
    EXPECT_EQUAL(std::string(::getenv("ANOTHER_VARIABLE")), "XYZ");
  }
}

class LocalEnvironment : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::LocalEnvironment";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_LOCALENVIRONMENT_H_
