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
#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "./TestConfig.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/ModelBiasCovariance.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ModBiasTestFixture : TestFixture {
 public:
  ModBiasTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "resolution");
    resol_.reset(new lorenz95::Resolution(res));
    biasconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ModelBias"));
    nobias_.reset(new eckit::LocalConfiguration());
    eckit::LocalConfiguration covarCfg(TestConfig::config(), "Covariance");
    lorenz95::ModelBiasCovariance covar(covarCfg, *resol_);
    bias1_ = biasconf_->getDouble("bias");
    bias2_ = 2.5 * bias1_;
    dbias_.reset(new lorenz95::ModelBiasCorrection(*resol_, covar.config()));
    dbias_->bias() = bias2_;
  }
  ~ModBiasTestFixture() {}
  boost::scoped_ptr<lorenz95::Resolution> resol_;
  boost::scoped_ptr<const eckit::LocalConfiguration> biasconf_;
  boost::scoped_ptr<const eckit::LocalConfiguration> nobias_;
  boost::scoped_ptr<lorenz95::ModelBiasCorrection> dbias_;
  double bias1_;
  double bias2_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_modelBias, ModBiasTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_constructor_bias) {
    lorenz95::ModelBias mbias(*resol_, *biasconf_);
    BOOST_CHECK_EQUAL(mbias.bias(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_constructor_nobias) {
    lorenz95::ModelBias mbias(*resol_, *nobias_);

    // because the biasconf_ is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    BOOST_CHECK_EQUAL(mbias.bias(), 0.0);
}
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_constructor_active) {
    lorenz95::ModelBias mbias1(*resol_, *biasconf_);
    lorenz95::ModelBias mbias2(*resol_, mbias1);
    BOOST_CHECK_EQUAL(mbias2.bias(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_constructor_inactive) {
    lorenz95::ModelBias mbias1(*resol_, *nobias_);
    lorenz95::ModelBias mbias2(*resol_, mbias1);

    // because the biasconf_ is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    BOOST_CHECK_EQUAL(mbias2.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_constructor_active_copy) {
    lorenz95::ModelBias mbias1(*resol_, *biasconf_);
    lorenz95::ModelBias mbias2(mbias1, true);
    BOOST_CHECK_EQUAL(mbias2.bias(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_constructor_active_nocopy) {
    lorenz95::ModelBias mbias1(*resol_, *biasconf_);
    lorenz95::ModelBias mbias2(mbias1, false);

    // because copy flag is set to false, bias is 0.0
    BOOST_CHECK_EQUAL(mbias2.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_constructor_inactive_copy) {
    lorenz95::ModelBias mbias1(*resol_, *nobias_);
    lorenz95::ModelBias mbias2(mbias1, true);

    // because the biasconf_ is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    BOOST_CHECK_EQUAL(mbias2.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_constructor_inactive_nocopy) {
    lorenz95::ModelBias mbias1(*resol_, *nobias_);
    lorenz95::ModelBias mbias2(mbias1, false);
    BOOST_CHECK_EQUAL(mbias2.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_destructor) {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_classname) {
    lorenz95::ModelBias mbias(*resol_, *biasconf_);
    BOOST_CHECK_EQUAL(mbias.classname(), "lorenz95::ModelBias");
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_compound_assignment_active) {
    lorenz95::ModelBias mbias(*resol_, *biasconf_);

    mbias += *dbias_;

    BOOST_CHECK_EQUAL(mbias.bias(), bias1_ + bias2_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_compound_assignment_inactive) {
    lorenz95::ModelBias mbias(*resol_, *nobias_);

    mbias += *dbias_;

    // because the biasconf_ is empty,
    // mbias bias value is not modified from 0.0
    BOOST_CHECK_EQUAL(mbias.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_bias) {
    lorenz95::ModelBias mbias(*resol_, *biasconf_);
    BOOST_CHECK_EQUAL(mbias.bias(), bias1_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_read) {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_write) {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBias_stream_output) {
    lorenz95::ModelBias mbias(*resol_, *biasconf_);

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ModelBiasTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << mbias;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double testBias = bias1_;
    double bias = 0.0;
    int biasStartPos = 12;  // length of "ModelBias = " is 12
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

BOOST_AUTO_TEST_SUITE_END()
}  // namespace test
