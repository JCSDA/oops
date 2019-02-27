/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <fstream>
#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>

#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/Resolution.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ResolutionTestFixture : TestFixtureBase<false> {
 public:
  ResolutionTestFixture() {
    testconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "resolution"));
  }
  ~ResolutionTestFixture() {}
  boost::scoped_ptr<const eckit::LocalConfiguration> testconf_;
};
// -----------------------------------------------------------------------------
CASE("test_resolution") {
  ResolutionTestFixture f;
// -----------------------------------------------------------------------------
  SECTION("test_resolution_constructor") {
    boost::scoped_ptr<lorenz95::Resolution> resol(new lorenz95::Resolution(*f.testconf_));
    EXPECT(resol.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_resolution_copy_constructor") {
    boost::scoped_ptr<lorenz95::Resolution> xx(new lorenz95::Resolution(*f.testconf_));
    boost::scoped_ptr<lorenz95::Resolution> resol(new lorenz95::Resolution(*xx));
    EXPECT(resol.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_resolution_get_npoints") {
    boost::scoped_ptr<lorenz95::Resolution> resol(new lorenz95::Resolution(*f.testconf_));
    EXPECT(resol->npoints() == f.testconf_->getInt("resol"));
  }
// -----------------------------------------------------------------------------
  SECTION("test_resolution_stream_output") {
    boost::scoped_ptr<lorenz95::Resolution> resol(new lorenz95::Resolution(*f.testconf_));

    // use the operator<< method to write the value to a file

    std::filebuf fb;
    std::string filename("ResolutionTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << *resol.get();
    fb.close();

    // then read the value that was written to the file
    int inputInt;
    std::string input;
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, input);

      try {
        inputInt = boost::lexical_cast<int>(input);
      }
      catch(boost::bad_lexical_cast const&) {
        // test fails because the value written to
        // the file can't be converted to an integer
        std::cout <<"operator<< incorrectly output a non-integer" << std::endl;
      }

      // it should equal the value that was used in the constructor
      EXPECT(inputInt == resol->npoints());
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      std::cout << "operator<< functionality cannot be determined" << std::endl;
    }
    inputFile.close();
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
