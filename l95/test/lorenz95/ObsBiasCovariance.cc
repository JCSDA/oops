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
#include <memory>


#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/ObsBiasCorrection.h"
#include "lorenz95/ObsBiasCovariance.h"
#include "lorenz95/Resolution.h"
#include "oops/runs/Test.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsBiasTestFixture : TestFixture {
 public:
  ObsBiasTestFixture() {
    biasconf_.reset(new eckit::LocalConfiguration(TestConfig::config(),
                                                  "Observations.Observation"));
    nobias_.reset(new eckit::LocalConfiguration());
  }
  ~ObsBiasTestFixture() {}
  std::unique_ptr<const eckit::LocalConfiguration> biasconf_;
  std::unique_ptr<const eckit::LocalConfiguration> nobias_;
};
// -----------------------------------------------------------------------------
CASE("test_obsBiasCovariance") {
  ObsBiasTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCovariance_constructor_conf") {
    lorenz95::ObsBiasCovariance obcovar(*fix.biasconf_);
    EXPECT(obcovar.active() == true);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCovariance_constructor_no_conf") {
    lorenz95::ObsBiasCovariance obcovar(*fix.nobias_);
    EXPECT(obcovar.active() == false);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCovariance_destructor") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCovariance_linearize") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCovariance_multiply_active") {
    lorenz95::ObsBiasCovariance obcovar(*fix.biasconf_);

    lorenz95::ObsBiasCorrection db1(*fix.biasconf_);
    db1.value() = 2.0;
    lorenz95::ObsBiasCorrection db2(db1, *fix.biasconf_);

    obcovar.multiply(db1, db2);

    const double stdev =
      fix.biasconf_->getDouble("ObsBiasCovariance.standard_deviation");
    EXPECT(db2.value() == db1.value() * stdev * stdev);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCovariance_multiply_inactive") {
    lorenz95::ObsBiasCovariance obcovar(*fix.nobias_);

    lorenz95::ObsBiasCorrection db1(*fix.nobias_);
    db1.value() = 2.0;
    lorenz95::ObsBiasCorrection db2(db1, *fix.biasconf_);

    obcovar.multiply(db1, db2);

    // because the OBC has empty config, the bias is set to 0.0
    EXPECT(db2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCovariance_invMult_active") {
    lorenz95::ObsBiasCovariance obcovar(*fix.biasconf_);

    lorenz95::ObsBiasCorrection db1(*fix.biasconf_);
    db1.value() = 2.0;
    lorenz95::ObsBiasCorrection db2(db1, *fix.biasconf_);

    obcovar.inverseMultiply(db1, db2);

    const double stdev =
      fix.biasconf_->getDouble("ObsBiasCovariance.standard_deviation");
    EXPECT(db2.value() == db1.value() * 1.0 / (stdev * stdev));
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCovariance_invMult_inactive") {
    lorenz95::ObsBiasCovariance obcovar(*fix.nobias_);

    lorenz95::ObsBiasCorrection db1(*fix.nobias_);
    db1.value() = 2.0;
    lorenz95::ObsBiasCorrection db2(db1, *fix.biasconf_);

    obcovar.inverseMultiply(db1, db2);

    // because the OBC has empty config, the bias is set to 0.0
    EXPECT(db2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
