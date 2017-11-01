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
#include "lorenz95/ObsBiasCorrection.h"
#include "lorenz95/ObsBiasCovariance.h"
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
BOOST_FIXTURE_TEST_SUITE(test_obsBiasCovariance, ObsBiasTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCovariance_constructor_conf) {
    lorenz95::ObsBiasCovariance obcovar(*covconf_);
    BOOST_CHECK_EQUAL(obcovar.active(), true);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCovariance_constructor_no_conf) {
    lorenz95::ObsBiasCovariance obcovar(*nobias_);
    BOOST_CHECK_EQUAL(obcovar.active(), false);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCovariance_destructor) {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCovariance_linearize) {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCovariance_multiply_active) {
    lorenz95::ObsBiasCovariance obcovar(*covconf_);

    lorenz95::ObsBiasCorrection db1(*covconf_);
    db1.value() = 2.0;
    lorenz95::ObsBiasCorrection db2(db1, *covconf_);

    obcovar.multiply(db1, db2);

    const double stdev = covconf_->getDouble("standard_deviation");
    BOOST_CHECK_EQUAL(db2.value(), db1.value() * stdev * stdev);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCovariance_multiply_inactive) {
    lorenz95::ObsBiasCovariance obcovar(*nobias_);

    lorenz95::ObsBiasCorrection db1(*nobias_);
    db1.value() = 2.0;
    lorenz95::ObsBiasCorrection db2(db1, *covconf_);

    obcovar.multiply(db1, db2);

    // because the OBC has empty config, the bias is set to 0.0
    BOOST_CHECK_EQUAL(db2.value(), 0.0);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCovariance_invMult_active) {
    lorenz95::ObsBiasCovariance obcovar(*covconf_);

    lorenz95::ObsBiasCorrection db1(*covconf_);
    db1.value() = 2.0;
    lorenz95::ObsBiasCorrection db2(db1, *covconf_);

    obcovar.inverseMultiply(db1, db2);

    const double stdev = covconf_->getDouble("standard_deviation");
    BOOST_CHECK_EQUAL(db2.value(), db1.value() * 1.0 / (stdev * stdev));
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_obsBiasCovariance_invMult_inactive) {
    lorenz95::ObsBiasCovariance obcovar(*nobias_);

    lorenz95::ObsBiasCorrection db1(*nobias_);
    db1.value() = 2.0;
    lorenz95::ObsBiasCorrection db2(db1, *covconf_);

    obcovar.inverseMultiply(db1, db2);

    // because the OBC has empty config, the bias is set to 0.0
    BOOST_CHECK_EQUAL(db2.value(), 0.0);
  }
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
}  // namespace test
