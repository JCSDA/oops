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
#include <boost/test/unit_test.hpp>

#include "./TestConfig.h"
#include "lorenz95/ObsBias.h"
#include "lorenz95/ObsBiasCorrection.h"
#include "eckit/config/LocalConfiguration.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsBiasTestFixture : TestFixture {
 public:
  ObsBiasTestFixture() {
    off_.reset(new eckit::LocalConfiguration());
    conf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ObsBiasCovariance"));
    eckit::LocalConfiguration bconf(TestConfig::config(), "ObsBias");
    obias_.reset(new lorenz95::ObsBias(bconf));
    bias1_ = bconf.getDouble("bias");
    bias2_ = 3.5 * bias1_;
    fact_ = 1.234;
  }
  ~ObsBiasTestFixture() {}
  boost::scoped_ptr<const eckit::LocalConfiguration> off_;
  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
  double bias1_;
  double bias2_;
  double fact_;
  boost::scoped_ptr<lorenz95::ObsBias> obias_;
};

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_obsBiasCorrection, ObsBiasTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_constructor_active) {
    lorenz95::ObsBiasCorrection dob(*conf_);
    BOOST_CHECK_EQUAL(dob.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_constructor_inactive) {
    lorenz95::ObsBiasCorrection dob(*off_);
    BOOST_CHECK_EQUAL(dob.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_copy_ctor_active_copy) {
    lorenz95::ObsBiasCorrection dob1(*conf_);
    dob1.value() = bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, true);

    BOOST_CHECK_EQUAL(dob2.value(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_copy_ctor_active_no_copy) {
    lorenz95::ObsBiasCorrection dob1(*conf_);
    dob1.value() = bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, false);

    // because the copy is false,
    // the active_ flag is true and the bias1_ value is 0.0
    BOOST_CHECK_EQUAL(dob2.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_copy_ctor_inactive_copy) {
    lorenz95::ObsBiasCorrection dob1(*off_);
    dob1.value() = bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, true);

    // because the cfg is empty when used,
    // the active_ flag is false and the bias1_ value is 0.0
    BOOST_CHECK_EQUAL(dob2.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_copy_ctor_inactive_no_copy) {
    lorenz95::ObsBiasCorrection dob1(*off_);
    dob1.value() = bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, false);

    // because the cfg is empty when used and the copy flag is false,
    // the active_ flag is false and the bias1_ value is 0.0
    BOOST_CHECK_EQUAL(dob2.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_copy_constructor_config) {
    lorenz95::ObsBiasCorrection dob1(*conf_);
    dob1.value() = bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, *conf_);

    BOOST_CHECK_EQUAL(dob2.value(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_copy_constructor_no_config) {
    lorenz95::ObsBiasCorrection dob1(*off_);
    dob1.value() = bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, eckit::LocalConfiguration());

    // because the covarCfg is empty when used (regardless of the cfg),
    // the active_ flag is false and the bias1_ value is 0.0
    BOOST_CHECK_EQUAL(dob2.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_diff_active) {
    lorenz95::ObsBias obias2(*conf_);
    obias2.value() = bias2_;

    lorenz95::ObsBiasCorrection dob(*conf_);

    dob.diff(*obias_, obias2);

    BOOST_CHECK_EQUAL(dob.value(), bias1_ - bias2_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_diff_inactive) {
    lorenz95::ObsBias obias2(*conf_);
    obias2.value() = bias2_;

    // construct the dob object with empty config
    lorenz95::ObsBiasCorrection dob(*off_);

    dob.diff(*obias_, obias2);

    // because the OBC has empty config the active_flag is false and
    // the diff will not be performed, leaving the OBC value unchanged
    BOOST_CHECK_EQUAL(dob.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_zero) {
    lorenz95::ObsBiasCorrection dob(*conf_);
    dob.value() = bias1_;

    dob.zero();

    BOOST_CHECK_EQUAL(dob.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_assignment_active) {
    lorenz95::ObsBiasCorrection dob1(*conf_);
    dob1.value() = bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = bias2_;

    dob1 = dob2;

    BOOST_CHECK_EQUAL(dob1.value(), bias2_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_assignment_inactive) {
    lorenz95::ObsBiasCorrection dob1(*off_);
    dob1.value() = bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = bias2_;

    dob1 = dob2;

    // because the OBC has empty config the active_ flag is false
    BOOST_CHECK_EQUAL(dob1.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_assignment_add_active) {
    lorenz95::ObsBiasCorrection dob1(*conf_);
    dob1.value() = bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = bias2_;

    dob1 += dob2;

    BOOST_CHECK_EQUAL(dob1.value(), bias1_ + bias2_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_assign_add_inactive) {
    lorenz95::ObsBiasCorrection dob1(*off_);
    dob1.value() = bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = bias2_;

    dob1 += dob2;

    // because the OBC has empty config, the bias value is unchanged
    BOOST_CHECK_EQUAL(dob1.value(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_assignment_subtract_active) {
    lorenz95::ObsBiasCorrection dob1(*conf_);
    dob1.value() = bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = bias2_;

    dob1 -= dob2;

    BOOST_CHECK_EQUAL(dob1.value(), bias1_ - bias2_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_assignment_subtract_inactive) {
    lorenz95::ObsBiasCorrection dob1(*off_);
    dob1.value() = bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = bias2_;

    dob1 -= dob2;

    // because the OBC has empty config, the bias value is unchanged
    BOOST_CHECK_EQUAL(dob1.value(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_assignment_multiply_active) {
    lorenz95::ObsBiasCorrection dob(*conf_);
    dob.value() = bias1_;

    dob *= fact_;

    BOOST_CHECK_EQUAL(dob.value(), bias1_ * fact_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_assignment_multiply_inactive) {
    lorenz95::ObsBiasCorrection dob(*off_);
    dob.value() = bias1_;

    dob *= fact_;

    // because the OBC has empty config, the bias value is unchanged
    BOOST_CHECK_EQUAL(dob.value(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_axpy_active) {
    lorenz95::ObsBiasCorrection dob1(*conf_);
    dob1.value() = bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = bias2_;

    dob1.axpy(fact_, dob2);

    BOOST_CHECK_EQUAL(dob1.value(), bias1_ + (fact_ * bias2_));
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_axpy_inactive) {
    lorenz95::ObsBiasCorrection dob1(*off_);
    dob1.value() = bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = bias2_;

    dob1.axpy(fact_, dob2);

    // because the OBC has empty config, the bias value is unchanged
    BOOST_CHECK_EQUAL(dob1.value(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_dot_product_with_active) {
    lorenz95::ObsBiasCorrection dob1(*conf_);
    dob1.value() = bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = bias2_;

    double dpwResult = dob1.dot_product_with(dob2);

    BOOST_CHECK_EQUAL(dpwResult, bias1_ * bias2_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_dot_product_with_inactive) {
    lorenz95::ObsBiasCorrection dob1(*off_);
    dob1.value() = bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = bias2_;

    double dpwResult = dob1.dot_product_with(dob2);

    // because the OBC has empty config, the result is 0
    BOOST_CHECK_EQUAL(dpwResult, 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_read) {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_write) {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCorrection_stream_output) {
    lorenz95::ObsBiasCorrection dob(*conf_);
    dob.value() = bias1_;

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ObsBiasCorrectionTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << dob;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double bias = 0.0;
    int biasStartPos = 20;  // length of "ObsBiasCorrection = " is 20
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, inputString);  // ignore first (blank) line
      getline(inputFile, inputString);

      inputBias = inputString.substr(biasStartPos);

      try {
        bias = boost::lexical_cast<double>(inputBias);
      }
      catch(boost::bad_lexical_cast const&) {
        BOOST_ERROR("operator<< incorrectly output a non-double");
      }

      BOOST_CHECK_CLOSE(bias1_, bias, 0.0001);
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
