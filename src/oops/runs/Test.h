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
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/testing/Test.h"
#include "eckit/utils/Tokenizer.h"

#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/CompareNVectors.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"

#include "test/TestEnvironment.h"

namespace oops {

// -----------------------------------------------------------------------------
class Test : public Application {
 public:
  explicit Test(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
  virtual ~Test() {}
  int execute(const eckit::Configuration & config, bool validate) const override;
 private:
  virtual void register_tests() const = 0;
  virtual std::string testid() const = 0;
  virtual void clear() const = 0;
  static bool init_unit_test() {return true;}
  std::string appname() const override {return "oops::Test running " + testid();}
};

// -----------------------------------------------------------------------------

int Test::execute(const eckit::Configuration & config, bool validate) const {
// Setup configuration for tests
  test::TestEnvironment::getInstance().setup(config);

// Generate the argc and argv arguments for unit_test_main(...)
  int argc = 1;
  char * argv[argc];
  char dummy[] = "abcde";
  argv[0] = dummy;

// Run the tests
  Log::trace() << "Registering the unit tests" << std::endl;
  register_tests();
  Log::trace() << "Running the unit tests" << std::endl;
  int result = eckit::testing::run_tests(argc, argv, false);
  Log::trace() << "Finished running the unit tests" << std::endl;
  Log::error() << "Finished running the unit tests, result = " << result << std::endl;

  this->clear();

// Return test status
  return result;
}

}  // namespace oops
#endif  // OOPS_RUNS_TEST_H_
