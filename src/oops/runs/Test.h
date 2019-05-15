/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_TEST_H_
#define OOPS_RUNS_TEST_H_

#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>
#include "eckit/config/Configuration.h"
#include "eckit/parser/Tokenizer.h"
#include "eckit/testing/Test.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"
#include "test/TestEnvironment.h"

namespace oops {

// -----------------------------------------------------------------------------
class Test : public Application {
 public:
  Test() {}
  virtual ~Test() {}
  int execute(const eckit::Configuration & config) const;
 private:
  virtual void register_tests() const = 0;
  virtual std::string testid() const = 0;
  static bool init_unit_test() {return true;}
  std::string appname() const {return "oops::Test running " + testid();}
};

// -----------------------------------------------------------------------------

int Test::execute(const eckit::Configuration & config) const {
// Setup configuration for tests
  test::TestEnvironment::getInstance().setup(config);

// Create a string version of argv
  std::vector<std::string> argvec;
  argvec.push_back(std::string("abcd"));

// Generate the argc and argv arguments for unit_test_main(...)
  int argc = argvec.size();
  char * argv[argc];
  for (int i = 0; i < argc; ++i) {
    argv[i] = new char[argvec[i].size()+1];
    strcpy(argv[i], argvec[i].c_str());
  }

// Run the tests
  Log::trace() << "Registering the unit tests" << std::endl;
  register_tests();
  Log::trace() << "Running the unit tests" << std::endl;
  int result = eckit::testing::run_tests(argc, argv, false);
  Log::trace() << "Finished running the unit tests" << std::endl;
  Log::error() << "Finished running the unit tests, result = " << result << std::endl;

// Tidy up
  for (int i = 0; i < argc; ++i) {
    delete [] argv[i];
  }

// Return test status
  return result;
}
// -----------------------------------------------------------------------------
template< typename T >
bool is_close(T a, T b, T epsilon) {
  // if nan or inf values, always return false
  if (std::isnan(a) || std::isnan(b) || std::isinf(a) || std::isinf(b)) return false;

  // otherwise, create a relative tolerance that is of the same type as a and b
  // (which is what is_approximately_equal wants) and call is_approximately_equal.
  T AbsA = fabs(a);
  T AbsB = fabs(b);
  T EpsAB = (AbsA < AbsB ? AbsB : AbsA) * epsilon;  // greater of AbsA, AbsB times epsilon
  bool test_status = eckit::types::is_approximately_equal(a , b, EpsAB);
  std::size_t num_digits = std::numeric_limits<T>::max_digits10;
  if (test_status) {
    Log::info() << "difference between " << std::setprecision(num_digits)
      << a << " and " << b << " is less than " << EpsAB << " (PASS)" << std::endl;
  } else {
    Log::info() << "difference between " << std::setprecision(num_digits)
      << a << " and " << b << " exceeds " << EpsAB << " (FAIL)" << std::endl;
  }

  return test_status;
}
  }  // namespace oops
#endif  // OOPS_RUNS_TEST_H_

