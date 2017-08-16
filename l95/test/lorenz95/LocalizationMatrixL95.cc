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
#include "lorenz95/StateL95.h"
#include "lorenz95/LocalizationMatrixL95.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class LocalizationTestFixture : TestFixture {
 public:
  LocalizationTestFixture() {
    file_.reset(new eckit::LocalConfiguration(TestConfig::config(), "state"));
    eckit::LocalConfiguration res(TestConfig::config(), "resolution");
    resol_.reset(new lorenz95::Resolution(res));
    date_str_ = "2014-09-12T09:35:00Z";
    time_.reset(new util::DateTime(date_str_));
    vars_.reset(new lorenz95::NoVariables(TestConfig::config()));
  }
  ~LocalizationTestFixture() {}
  boost::scoped_ptr<const eckit::LocalConfiguration> file_;
  boost::scoped_ptr<lorenz95::Resolution> resol_;
  std::string date_str_;
  boost::scoped_ptr<util::DateTime> time_;
  boost::scoped_ptr<lorenz95::NoVariables> vars_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_localizationMatrixL95,LocalizationTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_localizationMatrixL95_constructor) {

    eckit::LocalConfiguration resolCfg(TestConfig::config(), "resolution");
    lorenz95::Resolution resol(resolCfg);
    lorenz95::StateL95 state(*resol_, *vars_, *time_);
    eckit::LocalConfiguration cfg(TestConfig::config(), "Covariance");

    boost::scoped_ptr<lorenz95::LocalizationMatrixL95> locmat(
        new lorenz95::LocalizationMatrixL95(state, cfg));

    BOOST_CHECK(locmat.get() != NULL);
  }
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

}  // namespace test
