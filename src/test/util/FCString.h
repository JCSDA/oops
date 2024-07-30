/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */

#ifndef TEST_UTIL_FCSTRING_H_
#define TEST_UTIL_FCSTRING_H_

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"
#include "test/TestEnvironment.h"
#include "test/util/Fortran.h"

namespace test {

// -----------------------------------------------------------------------------
/*! Test Fixture for Fortran-C Interface Tools */

class FCStringFixture : private boost::noncopyable {
 public:
  static const eckit::Configuration & test() {return *getInstance().test_;}

 private:
  static FCStringFixture & getInstance() {
    static FCStringFixture theFCStringFixture;
    return theFCStringFixture;
  }

  FCStringFixture() {
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "test_f_c_string"));
  }

  ~FCStringFixture() {}

  std::unique_ptr<const eckit::LocalConfiguration> test_;
};

// -----------------------------------------------------------------------------
/*! Test transfer of string array from Fortran to C++
 *
 * \details This is intended to test the **f_c_push_string_vector()** utility
 * for passing a string vector from Fortran to C++
 *
 */

void testPushStringVector() {
  typedef FCStringFixture Test_;

  const eckit::Configuration * config = &Test_::test();

  std::vector<std::string> vec;
  test_push_string_vector_f(&config, vec);

  // retrieve directly from config file for comparison
  std::vector<std::string> vec_check(config->getStringVector("string_vec"));

  EXPECT(vec.size() == vec_check.size());
  for (std::size_t jj = 0; jj < vec_check.size(); ++jj) {
    EXPECT(vec[jj] == vec_check[jj]);
  }
}

// -----------------------------------------------------------------------------

class FCString : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~FCString() {}

 private:
  std::string testid() const override {return "test::FCString";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("util/FCString/testPushStringVector")
      { testPushStringVector(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------
}  // namespace test

#endif  // TEST_UTIL_FCSTRING_H_

