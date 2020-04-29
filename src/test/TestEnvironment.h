/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_TESTENVIRONMENT_H_
#define TEST_TESTENVIRONMENT_H_

#include <memory>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"

namespace test {

/// TestEnvironment is a singleton that defines the unit testing enviroment
/*! TestEnvironment contains globally available information for the unit
 *  tests. It is needed because there is no easy method to pass configuration
 *  and other data to the tests. By defining a singleton, each test can
 *  use TestEnvironment::getInstance() to get access to the global data.
 */

class TestEnvironment : private boost::noncopyable {
 public:
  static TestEnvironment & getInstance() {
    static TestEnvironment theTestEnvironment;
    return theTestEnvironment;
  }

  void setup(const eckit::Configuration & conf) {
    config_.reset(new eckit::LocalConfiguration(conf));
  }

  static const eckit::Configuration & config() {return *getInstance().config_;}

 private:
  TestEnvironment() {}
  ~TestEnvironment() {}

  std::unique_ptr<const eckit::Configuration> config_;
};

}  // namespace test

#endif  // TEST_TESTENVIRONMENT_H_
