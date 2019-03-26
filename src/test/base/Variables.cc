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
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Variables.h"
#include "test/TestEnvironment.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class VariablesFixture : TestFixture {
 public:
  VariablesFixture() {
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "Variables"));
  }

  ~VariablesFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
};
// -----------------------------------------------------------------------------
CASE("VariablesFixture") {
  VariableFixture fix;
// -----------------------------------------------------------------------------

SECTION("testConstructor") {
  boost::scoped_ptr<oops::Variables> vars(new oops::Variables(*fix.conf_));
  EXPECT(vars.get());

  vars.reset();
  EXPECT(!vars.get());
}

// -----------------------------------------------------------------------------

SECTION("testCopyConstructor") {
  boost::scoped_ptr<oops::Variables> vars(new oops::Variables(*fix.conf_));

  boost::scoped_ptr<oops::Variables> other(new oops::Variables(*vars));
  EXPECT(other.get());

  other.reset();
  EXPECT(!other.get());

  EXPECT(vars.get());
}

// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------

}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
