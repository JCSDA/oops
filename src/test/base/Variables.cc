/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <string>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
class VariablesFixture {
 public:
  VariablesFixture() {
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "Variables"));
  }

  ~VariablesFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_Variables, VariablesFixture)
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testConstructor) {
  boost::scoped_ptr<oops::Variables> vars(new oops::Variables(*conf_));
  BOOST_CHECK(vars.get());

  vars.reset();
  BOOST_CHECK(!vars.get());
}

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testCopyConstructor) {
  boost::scoped_ptr<oops::Variables> vars(new oops::Variables(*conf_));

  boost::scoped_ptr<oops::Variables> other(new oops::Variables(*vars));
  BOOST_CHECK(other.get());

  other.reset();
  BOOST_CHECK(!other.get());

  BOOST_CHECK(vars.get());
}

// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------

}  // namespace test
