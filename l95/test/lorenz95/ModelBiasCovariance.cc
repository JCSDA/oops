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

#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "./TestConfig.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/ModelBiasCovariance.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ModBiasCovTestFixture : TestFixture {
 public:
  ModBiasCovTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "resolution");
    resol_.reset(new lorenz95::Resolution(res));
    covconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ModelBiasCovariance"));
    nobias_.reset(new eckit::LocalConfiguration());
  }
  ~ModBiasCovTestFixture() {}
  boost::scoped_ptr<lorenz95::Resolution> resol_;
  boost::scoped_ptr<const eckit::LocalConfiguration> covconf_;
  boost::scoped_ptr<const eckit::LocalConfiguration> nobias_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_modelBiasCovariance, ModBiasCovTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCovariance_constructor_conf) {
    lorenz95::ModelBiasCovariance bcovar(*covconf_, *resol_);
    BOOST_CHECK_EQUAL(bcovar.active(), true);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCovariance_constructor_no_conf) {
    lorenz95::ModelBiasCovariance bcovar(*nobias_, *resol_);
    BOOST_CHECK_EQUAL(bcovar.active(), false);
}
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCovariance_linearize) {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCovariance_multiply_active) {
    // construct the ModelBiasCorrection object
    lorenz95::ModelBiasCovariance bcovar(*covconf_, *resol_);
    lorenz95::ModelBiasCorrection dbias1(*resol_, *covconf_);
    dbias1.bias() = 2.0;
    lorenz95::ModelBiasCorrection dbias2(dbias1, true);

    bcovar.multiply(dbias1, dbias2);

    double stdev = covconf_->getDouble("standard_deviation");
    BOOST_CHECK_EQUAL(dbias2.bias(), dbias1.bias() * stdev * stdev);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCovariance_multiply_inactive) {
    // construct the ModelBiasCorrection object
    lorenz95::ModelBiasCovariance bcovar(*nobias_, *resol_);
    lorenz95::ModelBiasCorrection dbias1(*resol_, *covconf_);
    dbias1.bias() = 2.0;
    lorenz95::ModelBiasCorrection dbias2(dbias1, true);

    bcovar.multiply(dbias1, dbias2);

    // because the covconf_ is empty, the bias is set to 0
    BOOST_CHECK_EQUAL(dbias2.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCovariance_invMult_active) {
    // construct the ModelBiasCorrection object
    lorenz95::ModelBiasCovariance bcovar(*covconf_, *resol_);
    lorenz95::ModelBiasCorrection dbias1(*resol_, *covconf_);
    dbias1.bias() = 2.0;
    lorenz95::ModelBiasCorrection dbias2(dbias1, true);

    bcovar.inverseMultiply(dbias1, dbias2);

    double stdev = covconf_->getDouble("standard_deviation");
    BOOST_CHECK_EQUAL(dbias2.bias(), dbias1.bias() *1.0 / (stdev * stdev));
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCovariance_invMult_inactive) {
    // construct the ModelBiasCorrection object
    lorenz95::ModelBiasCovariance bcovar(*nobias_, *resol_);
    lorenz95::ModelBiasCorrection dbias1(*resol_, *covconf_);
    dbias1.bias() = 2.0;
    lorenz95::ModelBiasCorrection dbias2(dbias1, true);

    bcovar.inverseMultiply(dbias1, dbias2);

    // because the covconf_ is empty, the bias is set to 0
    BOOST_CHECK_EQUAL(dbias2.bias(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_modelBiasCovariance_active) {
    lorenz95::ModelBiasCovariance bcovar(*covconf_, *resol_);
    BOOST_CHECK_EQUAL(bcovar.active(), true);
  }
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
}  // namespace test
