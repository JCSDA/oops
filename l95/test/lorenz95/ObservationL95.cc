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
#include "lorenz95/ObservationL95.h"
#include "lorenz95/ObsTable.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsTestFixture : TestFixture {
 public:
  ObsTestFixture() {
    const eckit::LocalConfiguration conf(TestConfig::config(), "Observations");
    const util::DateTime bgn(conf.getString("window_begin"));
    const util::DateTime end(conf.getString("window_end"));
    const eckit::LocalConfiguration otconf(conf, "Observation");
    ot_.reset(new lorenz95::ObsTable(otconf, bgn, end));
  }
  ~ObsTestFixture() {}
  boost::scoped_ptr<lorenz95::ObsTable> ot_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_ObsL95, ObsTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsL95_constructor) {
    boost::scoped_ptr<lorenz95::ObservationL95> obs(lorenz95::ObservationL95::create(*ot_, TestConfig::config()));
    BOOST_CHECK(obs.get() != NULL);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_observationL95_classname) {
    boost::scoped_ptr<lorenz95::ObservationL95> obs(lorenz95::ObservationL95::create(*ot_, TestConfig::config()));
    BOOST_CHECK_EQUAL(obs->classname(), "lorenz95::ObservationL95");
  }
// -----------------------------------------------------------------------------
//   BOOST_AUTO_TEST_CASE(test_observationL95_obsEquiv) {
//     // create an ObservationL95 object
//     boost::shared_ptr<const eckit::LocalConfiguration> fullFileCfg(
//         new eckit::LocalConfiguration("test.xml", new util::XmlDom()));
//     boost::scoped_ptr<const eckit::LocalConfiguration> cfg(
//         new eckit::LocalConfiguration(TestConfig::config(),
//                          "Observations/Observation",
//                          true));
//     boost::scoped_ptr<lorenz95::ObservationL95>
//         observationL95(lorenz95::ObservationL95::create(*cfg));
//
//     // create a gomL95 object
//     int vecSize = 5;
//     std::vector<int> intVec;
//     for(int i = 0; i < vecSize; ++i) {
//       intVec.push_back(i);
//     }
//     boost::scoped_ptr<lorenz95::GomL95> gomL95(new lorenz95::GomL95(intVec));
//     // populate the locval_ vector with values from 1.1 to 5.5
//     for(int i = 0; i < vecSize; ++i){
//       (*gomL95)[i] = ((i + 1) * 1.1);
//     }
//
//     // create an ObsVec1D object
//     boost::scoped_ptr<lorenz95::ObsVec1D> obsVec1D(
//         new lorenz95::ObsVec1D(vecSize));
//
//     // create an ObsBias object
//     boost::scoped_ptr<const eckit::LocalConfiguration> biasCfg(
//         new eckit::LocalConfiguration(TestConfig::config(), "ObsBias", true));
//     boost::scoped_ptr<lorenz95::ObsBias> obsBias(
//         new lorenz95::ObsBias(*biasCfg));
//
//     observationL95->obsEquiv(*gomL95, *obsVec1D, *obsBias);
//
//     for(int i = 0; i < gomL95->nobs(); ++i) {
//       BOOST_CHECK_EQUAL((*obsVec1D)(gomL95->getindx(i)), (*gomL95)[i] + obsBias->value());
//     }
//   }
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------

}  // namespace test
