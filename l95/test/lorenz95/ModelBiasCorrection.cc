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

#include "eckit/config/LocalConfiguration.h"
#include "./TestConfig.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ModBiasTestFixture : TestFixture {
 public:
  ModBiasTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "resolution");
    resol_.reset(new lorenz95::Resolution(res));
    conf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ModelBiasCovariance"));
    nobias_.reset(new eckit::LocalConfiguration());
    bias1_ = TestConfig::config().getDouble("ModelBias.bias");
    bias2_ = 2.5 * bias1_;
    fact_ = 1.2345;
  }
  ~ModBiasTestFixture() {}
  boost::scoped_ptr<lorenz95::Resolution> resol_;
  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
  boost::scoped_ptr<const eckit::LocalConfiguration> nobias_;
  double bias1_;
  double bias2_;
  double fact_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_modelBiasCorrection, ModBiasTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_constructor_config) {
    boost::scoped_ptr<lorenz95::ModelBiasCorrection> dx(
      new lorenz95::ModelBiasCorrection(*resol_, *conf_));
    BOOST_CHECK(dx.get() != NULL);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_constructor_no_config) {
    boost::scoped_ptr<lorenz95::ModelBiasCorrection> dx(
      new lorenz95::ModelBiasCorrection(*resol_, *nobias_));
    BOOST_CHECK(dx.get() != NULL);
}
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_copy_ctor_active_copy) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *conf_);
    dx1.bias() = bias1_;

    lorenz95::ModelBiasCorrection dx2(dx1, true);

    BOOST_CHECK_EQUAL(dx2.bias(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_copy_ctor_active_no_copy) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *conf_);
    dx1.bias() = bias1_;

    // construct a copy of it with the copy flag set to false
    lorenz95::ModelBiasCorrection dx2(dx1, false);

    // because the copy is false,
    // the active_ flag is true and the bias_ value is 0.0
    BOOST_CHECK_EQUAL(dx2.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_copy_ctor_inactive_copy) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *nobias_);
    dx1.bias() = bias1_;

    // construct a copy of it with the copy flag set to true
    lorenz95::ModelBiasCorrection dx2(dx1, true);

    // because the cfg is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    BOOST_CHECK_EQUAL(dx2.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_copy_ctor_inactive_no_copy) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *nobias_);
    dx1.bias() = bias1_;

    // construct a copy of it with the copy flag set to false
    lorenz95::ModelBiasCorrection dx2(dx1, false);

    // because the cfg is empty when used and the copy flag is false,
    // the active_ flag is false and the bias_ value is 0.0
    BOOST_CHECK_EQUAL(dx2.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_copy_ctor_config_active) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *conf_);
    dx1.bias() = bias1_;

    lorenz95::ModelBiasCorrection dx2(dx1, *conf_);

    BOOST_CHECK_EQUAL(dx2.bias(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_copy_ctor_config_inactive) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *nobias_);
    dx1.bias() = bias1_;

    lorenz95::ModelBiasCorrection dx2(dx1, *conf_);

     // because the covarCfg is empty when used (regardless of the cfg),
     // the active_ flag is false and the bias_ value is 0.0
     BOOST_CHECK_EQUAL(dx2.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_diff_active) {
    lorenz95::ModelBias xx1(*resol_, *conf_);
    xx1.bias() = bias1_;
    lorenz95::ModelBias xx2(*resol_, *conf_);
    xx2.bias() = bias2_;

    lorenz95::ModelBiasCorrection dx(*resol_, *conf_);

    dx.diff(xx1, xx2);

    BOOST_CHECK_EQUAL(dx.bias(), bias1_ - bias2_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_diff_inactive) {
    lorenz95::ModelBias xx1(*resol_, *conf_);
    xx1.bias() = bias1_;
    lorenz95::ModelBias xx2(*resol_, *conf_);
    xx2.bias() = bias2_;

    lorenz95::ModelBiasCorrection dx(*resol_, *nobias_);

    dx.diff(xx1, xx2);

    // because the active_ flag is false, the bias cannot be updated
    BOOST_CHECK_EQUAL(dx.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_zero) {
    lorenz95::ModelBiasCorrection dx(*resol_, *conf_);
    dx.bias() = bias1_;

    dx.zero();

    BOOST_CHECK_EQUAL(dx.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_assignment_active) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *conf_);
    dx1.bias() = bias1_;
    lorenz95::ModelBiasCorrection dx2(*resol_, *conf_);
    dx2.bias() = bias2_;

    dx1 = dx2;

    // the original MBC should have the same bias value as the copy MBC
    BOOST_CHECK_EQUAL(dx1.bias(), bias2_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_assignment_inactive) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *nobias_);
    dx1.bias() = bias1_;

    lorenz95::ModelBiasCorrection dx2(*resol_, *conf_);
    dx2.bias() = bias2_;

    dx1 = dx2;

    // the active_ value is zero, so the bias will be zero
    BOOST_CHECK_EQUAL(dx1.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_assignment_add_active) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *conf_);
    dx1.bias() = bias1_;
    lorenz95::ModelBiasCorrection dx2(*resol_, *conf_);
    dx2.bias() = bias2_;

    dx1 += dx2;

    BOOST_CHECK_EQUAL(dx1.bias(), bias1_ + bias2_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_assignment_add_inactive) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *nobias_);
    dx1.bias() = bias1_;
    lorenz95::ModelBiasCorrection dx2(*resol_, *nobias_);
    dx2.bias() = bias2_;

    dx1 += dx2;

    // the active_ value is zero, so the bias will be unchanged
    BOOST_CHECK_EQUAL(dx1.bias(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_assignment_subtract_active) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *conf_);
    dx1.bias() = bias1_;
    lorenz95::ModelBiasCorrection dx2(*resol_, *conf_);
    dx2.bias() = bias2_;

    dx1 -= dx2;

    BOOST_CHECK_EQUAL(dx1.bias(), bias1_ - bias2_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_assignment_subtract_inactive) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *nobias_);
    dx1.bias() = bias1_;
    lorenz95::ModelBiasCorrection dx2(*resol_, *nobias_);
    dx2.bias() = bias2_;

    dx1 -= dx2;

    // the active_ value is zero, so the bias will be unchanged
    BOOST_CHECK_EQUAL(dx1.bias(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_assignment_multiply_active) {
    lorenz95::ModelBiasCorrection dx(*resol_, *conf_);
    dx.bias() = bias1_;

    dx *= fact_;

    BOOST_CHECK_EQUAL(dx.bias(), bias1_ * fact_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_assignment_multiply_inactive) {
    lorenz95::ModelBiasCorrection dx(*resol_, *nobias_);
    dx.bias() = bias1_;

    dx *= fact_;

    // the active_ value is zero, so the bias will be unchanged
    BOOST_CHECK_EQUAL(dx.bias(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_axpy_active) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *conf_);
    dx1.bias() = bias1_;
    lorenz95::ModelBiasCorrection dx2(*resol_, *conf_);
    dx2.bias() = bias2_;

    dx1.axpy(fact_, dx2);

    BOOST_CHECK_EQUAL(dx1.bias(), (bias1_ + fact_ * bias2_));
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_axpy_inactive) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *nobias_);
    dx1.bias() = bias1_;
    lorenz95::ModelBiasCorrection dx2(*resol_, *nobias_);
    dx2.bias() = bias2_;

    dx1.axpy(fact_, dx2);

    // the active_ value is zero, so the bias will be unchanged
    BOOST_CHECK_EQUAL(dx1.bias(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_dot_product_with_active) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *conf_);
    dx1.bias() = bias1_;
    lorenz95::ModelBiasCorrection dx2(*resol_, *conf_);
    dx2.bias() = bias2_;

    double dpwResult = dx1.dot_product_with(dx2);

    BOOST_CHECK_EQUAL(dpwResult, bias1_ * bias2_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_dot_product_with_inactive) {
    lorenz95::ModelBiasCorrection dx1(*resol_, *nobias_);
    dx1.bias() = bias1_;
    lorenz95::ModelBiasCorrection dx2(*resol_, *nobias_);
    dx2.bias() = bias2_;

    double dpwResult = dx1.dot_product_with(dx2);

    // because of the empty config, the result is 0.0
    BOOST_CHECK_EQUAL(dpwResult, 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_read) {
    // nothing to test
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_write) {
    // nothing to test
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_stream_output) {
    lorenz95::ModelBiasCorrection dx(*resol_, *conf_);
    dx.bias() = bias1_;

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ModelBiasCorrectionTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << dx;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double testBias = bias1_;
    double bias = 0.0;
    int biasStartPos = 22;  // length of "ModelBiasCorrection = " is 22
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

      BOOST_CHECK_CLOSE(testBias, bias, 0.0001);
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      BOOST_ERROR("operator<< functionality cannot be determined");
    }
    inputFile.close();
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCorrection_bias) {
    lorenz95::ModelBiasCorrection dx(*resol_, *conf_);
    dx.bias() = bias1_;

    // this one test checks both the setting and getting of the bias value
    BOOST_CHECK_EQUAL(dx.bias(), bias1_);
  }
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
}  // namespace test
