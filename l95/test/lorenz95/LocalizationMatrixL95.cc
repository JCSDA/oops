/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "./TestConfig.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/LocalizationMatrixL95.h"
#include "test/TestFixture.h"

namespace test {

BOOST_FIXTURE_TEST_SUITE(test_localizationMatrixL95, TestFixture)

// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_localizationMatrixL95_constructor) {
    eckit::LocalConfiguration resolCfg(TestConfig::config(), "resolution");
    lorenz95::Resolution resol(resolCfg);
    eckit::LocalConfiguration cfg(TestConfig::config(), "Covariance");

    boost::scoped_ptr<lorenz95::LocalizationMatrixL95> locmat(
        new lorenz95::LocalizationMatrixL95(resol, cfg));

    BOOST_CHECK(locmat.get() != NULL);
  }
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

}  // namespace test
