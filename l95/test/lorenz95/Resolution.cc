/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "./TestConfig.h"
#include "lorenz95/Resolution.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ResolutionTestFixture : TestFixture {
 public:
  ResolutionTestFixture() {
    testconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "resolution"));
  }
  ~ResolutionTestFixture() {}
  boost::scoped_ptr<const eckit::LocalConfiguration> testconf_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_resolution, ResolutionTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_resolution_constructor) {
    boost::scoped_ptr<lorenz95::Resolution> resol(new lorenz95::Resolution(*testconf_));
    BOOST_CHECK(resol.get() != NULL);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_resolution_copy_constructor) {
    boost::scoped_ptr<lorenz95::Resolution> xx(new lorenz95::Resolution(*testconf_));
    boost::scoped_ptr<lorenz95::Resolution> resol(new lorenz95::Resolution(*xx));
    BOOST_CHECK(resol.get() != NULL);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_resolution_get_npoints) {
    boost::scoped_ptr<lorenz95::Resolution> resol(new lorenz95::Resolution(*testconf_));
    BOOST_CHECK_EQUAL(resol->npoints(), testconf_->getInt("resol"));
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_resolution_stream_output) {
    boost::scoped_ptr<lorenz95::Resolution> resol(new lorenz95::Resolution(*testconf_));

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
        BOOST_ERROR("operator<< incorrectly output a non-integer");
      }

      // it should equal the value that was used in the constructor
      BOOST_CHECK_EQUAL(inputInt, resol->npoints());
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      BOOST_ERROR("operator<< functionality cannot be determined");
    }
    inputFile.close();
  }
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------

}  // namespace test
