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

#include <boost/test/unit_test.hpp>

#include "oops/base/Variables.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"
#include "eckit/config/Configuration.h"

namespace {

// -----------------------------------------------------------------------------
class VariablesFixture : private boost::noncopyable {
 public:
  static const eckit::Configuration & getConfig() {return *getInstance().conf_;}

 private:
  static VariablesFixture & getInstance() {
    static VariablesFixture theVariablesFixture;
    return theVariablesFixture;
  }

  VariablesFixture() {
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "Variables"));
  }

  ~VariablesFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
};

// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(test_Variables)
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testConstructor) {
  boost::scoped_ptr<Variables> vars(new Variables(VariablesFixture::getConfig()));
  BOOST_CHECK(vars.get());

  vars.reset();
  BOOST_CHECK(!vars.get());
}

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testCopyConstructor) {
  boost::scoped_ptr<Variables> vars(new Variables(VariablesFixture::getConfig()));

  boost::scoped_ptr<Variables> other(new Variables(*vars));
  BOOST_CHECK(other.get());

  other.reset();
  BOOST_CHECK(!other.get());

  BOOST_CHECK(vars.get());
}

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
