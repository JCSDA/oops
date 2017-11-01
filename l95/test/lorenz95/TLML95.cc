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

#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "./TestConfig.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/ModelBiasCovariance.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"
#include "lorenz95/TLML95.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class TlmTestFixture : TestFixture {
 public:
  TlmTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "resolution");
    resol_.reset(new lorenz95::Resolution(res));
    eckit::LocalConfiguration mod(TestConfig::config(), "model");
    model_.reset(new lorenz95::ModelL95(*resol_, mod));
    tlconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "linearmodel"));
  }
  ~TlmTestFixture() {}
  boost::scoped_ptr<lorenz95::ModelL95>   model_;
  boost::scoped_ptr<lorenz95::Resolution> resol_;
  boost::scoped_ptr<const eckit::LocalConfiguration> tlconf_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_tlmL95, TlmTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_tlmL95_constructor) {
    boost::scoped_ptr<lorenz95::TLML95> tlm(new lorenz95::TLML95(*resol_, *tlconf_));
    BOOST_CHECK(tlm.get() != NULL);
  }
// -----------------------------------------------------------------------------
/*
  BOOST_AUTO_TEST_CASE(test_tlmL95_set_get_Trajectory) {
    lorenz95::TLML95 tlm(*resol_, *tlconf_, *model_);

    // construct a StateL95 object
    std::string date_string("2014-09-12T09:35:00Z");
    util::DateTime dt(date_string);
    oops::Variables vars(util::emptyConfig());
    boost::scoped_ptr<lorenz95::StateL95> stateL95(new lorenz95::StateL95(*resol_, vars, dt));

    // construct a ModelBias object
    boost::shared_ptr<const eckit::LocalConfiguration> biasCfg(new eckit::LocalConfiguration(TestConfig::config(), "ModelBias"));
    boost::scoped_ptr<lorenz95::ModelBias> modelBias(new lorenz95::ModelBias(*biasCfg, *resol_));

    tlm.setTrajectory(*stateL95, *modelBias);

    boost::scoped_ptr<lorenz95::ModelTrajectory> modelTraj(new lorenz95::ModelTrajectory(true));
    std::cout << "PMC: getTraj result <" << tlm.getTrajectory(dt)->get(1) << ">" << std::endl;
    modelTraj->set(tlm.getTrajectory(dt)->get(1));
    //std::cout << "PMC: modelTraj 1st " << modelTraj->get(1)
    std::cout << "PMC: modelTraj resol  <" << modelTraj->get(1).resol() << ">" << std::endl;
    for(int i = 0; i < modelTraj->get(1).resol(); ++i) {
      std::cout << "PMC: modelTraj FieldL95 elements " << modelTraj->get(1)[i] << std::endl;
    }
  }
*/
// -----------------------------------------------------------------------------
/*
  BOOST_AUTO_TEST_CASE(test_tlmL95_stepTL) {
    lorenz95::TLML95 tlm(*resol_, *tlconf_, *model_);

    // construct a FieldL95 object
    boost::scoped_ptr<lorenz95::FieldL95> fieldL95(new lorenz95::FieldL95(*resol_));

    // construct a datetime object
    std::string dateString("2014-10-08T11:25:45Z");
    util::DateTime dt(dateString);

    //construct a ModelBiasCorrection object
    boost::scoped_ptr<const eckit::LocalConfiguration> covarCfg(new eckit::LocalConfiguration(TestConfig::config(), "Covariance"));
    boost::scoped_ptr<lorenz95::ModelBiasCovariance> modelBiasCovariance(
        new lorenz95::ModelBiasCovariance(*covarCfg, *resol_));
    boost::scoped_ptr<lorenz95::ModelBiasCorrection> modelBiasCorrection(
        new lorenz95::ModelBiasCorrection(*modelBiasCovariance));

    tlm.stepTL(*fieldL95, dt, *modelBiasCorrection);

    std::cout << "PMC: original dt was " << dateString << std::endl;
    std::cout << "PMC: modified dt is  " << dt.toString() << std::endl;
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_tlmL95_stepAD) {
  }
*/
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_tlmL95_get_classname) {
    lorenz95::TLML95 tlm(*resol_, *tlconf_);
    BOOST_CHECK_EQUAL(tlm.classname(), "lorenz95::TLML95");
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_tlmL95_get_timestep) {
    lorenz95::TLML95 tlm(*resol_, *tlconf_);
    util::Duration dt(tlconf_->getString("tstep"));
    BOOST_CHECK_EQUAL(tlm.timeResolution().toSeconds(), dt.toSeconds());
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_tlmL95_get_resolution) {
    eckit::LocalConfiguration rescf(TestConfig::config(), "resolution");
    lorenz95::Resolution resol(rescf);
    lorenz95::TLML95 tlm(*resol_, *tlconf_);
    BOOST_CHECK_EQUAL(tlm.resolution().npoints(), resol.npoints());
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_tlmL95_stream_output) {
    lorenz95::TLML95 tlm(*resol_, *tlconf_);

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("TLML95Test.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << tlm;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputStartsWith;
    std::string testString("L95 Model Trajectory, nstep=");
    int endPos = testString.size();
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, inputString);

      inputStartsWith = inputString.substr(0, endPos);

      BOOST_CHECK_EQUAL(inputStartsWith, testString);
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
