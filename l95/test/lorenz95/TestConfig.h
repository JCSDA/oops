/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef L95_TEST_LORENZ95_TESTCONFIG_H_
#define L95_TEST_LORENZ95_TESTCONFIG_H_

#include <memory>
#include <string>


#include "eckit/config/LocalConfiguration.h"
#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/runtime/Main.h"

namespace test {

// -----------------------------------------------------------------------------

class TestConfig {
 public:
  static const eckit::Configuration & config() {
    static TestConfig test;
    return *test.config_;
  }

  TestConfig() {
    ASSERT(eckit::Main::ready());

    int narg = eckit::Main::instance().argc();
    ASSERT(narg >= 2);
    eckit::PathName fname = eckit::Main::instance().argv(narg-1);
    config_.reset(new eckit::YAMLConfiguration(fname));
  }

  ~TestConfig() {}

 private:
  std::unique_ptr<const eckit::YAMLConfiguration> config_;
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // L95_TEST_LORENZ95_TESTCONFIG_H_
