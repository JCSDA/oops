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
#include "lorenz95/ObsBias.h"
#include "lorenz95/ObsBiasCorrection.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsBiasTestFixture : TestFixture {
 public:
  ObsBiasTestFixture() {
    biasconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ObsBias"));
    nobias_.reset(new eckit::LocalConfiguration());
    covconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "Covariance"));
  }
  ~ObsBiasTestFixture() {}
  boost::scoped_ptr<const eckit::LocalConfiguration> biasconf_;
  boost::scoped_ptr<const eckit::LocalConfiguration> nobias_;
  boost::scoped_ptr<const eckit::LocalConfiguration> covconf_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_obsBias, ObsBiasTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_constructor_bias) {
    lorenz95::ObsBias ob(*biasconf_);
    BOOST_CHECK_EQUAL(ob.value(), biasconf_->getDouble("bias"));
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_constructor_no_bias) {
    lorenz95::ObsBias ob(*nobias_);
    BOOST_CHECK_EQUAL(ob.value(), 0.0);
}
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_copy_constructor_config_copy) {
    lorenz95::ObsBias ob1(*biasconf_);
    lorenz95::ObsBias ob2(ob1, true);
    BOOST_CHECK_EQUAL(ob2.value(), biasconf_->getDouble("bias"));
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_copy_constructor_config_no_copy) {
    lorenz95::ObsBias ob1(*biasconf_);
    lorenz95::ObsBias ob2(ob1, false);

    // bias value is 0 because copy flag was set to false
    BOOST_CHECK_EQUAL(ob2.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_copy_constructor_no_config_copy) {
    lorenz95::ObsBias ob1(*nobias_);
    lorenz95::ObsBias ob2(ob1, true);

    // bias value is 0 because an empty config was used
    BOOST_CHECK_EQUAL(ob2.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_copy_constructor_no_config_no_copy) {
    lorenz95::ObsBias ob1(*nobias_);
    lorenz95::ObsBias ob2(ob1, true);

    // bias value is 0 because an empty config was used
    BOOST_CHECK_EQUAL(ob2.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_destructor) {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_compound_assignment_add_active) {
    lorenz95::ObsBias ob(*biasconf_);

    // construct an obsBiasCorrection object
    lorenz95::ObsBiasCorrection obsBiasCorrection(*covconf_);
    obsBiasCorrection.value() = 3.14;

    ob += obsBiasCorrection;

    BOOST_CHECK_EQUAL(ob.value(), biasconf_->getDouble("bias") + 3.14);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_compound_assignment_add_inactive) {
    lorenz95::ObsBias ob(*nobias_);

    // construct an obsBiasCorrection object
    lorenz95::ObsBiasCorrection obsBiasCorrection(*covconf_);
    obsBiasCorrection.value() = 3.14;

    ob += obsBiasCorrection;

    // the bias value will be 0 because the ob had an empty config
    BOOST_CHECK_EQUAL(ob.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_value) {
    lorenz95::ObsBias ob(*biasconf_);
    BOOST_CHECK_EQUAL(ob.value(), biasconf_->getDouble("bias"));
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_read) {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_write) {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBias_stream_output) {
    lorenz95::ObsBias ob(*biasconf_);

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ObsBiasTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << ob;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double testBias = biasconf_->getDouble("bias");
    double bias = 0.0;
    int biasStartPos = 10;  // length of "ObsBias = " is 10
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, inputString);
      getline(inputFile, inputString);

      inputBias = inputString.substr(biasStartPos);

      try {
        bias = boost::lexical_cast<double>(inputBias);
      }
      catch(boost::bad_lexical_cast const&) {
        BOOST_ERROR("operator<< incorrectly output a non-double");
      }

      BOOST_CHECK_CLOSE(testBias, bias, 0.0001);
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      BOOST_ERROR("operator<< functionality cannot be determined");
    }
    inputFile.close();
  }
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
}  // namespace test
