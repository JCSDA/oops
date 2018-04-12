/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "./TestConfig.h"
#include "lorenz95/FieldL95.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"
#include "util/DateTime.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class StateTestFixture : TestFixture {
 public:
  StateTestFixture() {
    file_.reset(new eckit::LocalConfiguration(TestConfig::config(), "state"));
    eckit::LocalConfiguration res(TestConfig::config(), "resolution");
    resol_.reset(new lorenz95::Resolution(res));
    date_str_ = file_->getString("date");
    time_.reset(new util::DateTime(date_str_));
    vars_.reset(new oops::Variables(TestConfig::config()));
  }
  ~StateTestFixture() {}
  boost::scoped_ptr<const eckit::LocalConfiguration> file_;
  boost::scoped_ptr<lorenz95::Resolution> resol_;
  std::string date_str_;
  boost::scoped_ptr<util::DateTime> time_;
  boost::scoped_ptr<oops::Variables> vars_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_StateL95, StateTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_stateL95_constructor) {
    boost::scoped_ptr<lorenz95::StateL95> xx(new lorenz95::StateL95(*resol_, *vars_, *time_));
    BOOST_CHECK(xx.get() != NULL);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_stateL95_readin_constructor) {
    boost::scoped_ptr<lorenz95::StateL95> xx(new lorenz95::StateL95(*resol_, *file_));
    BOOST_CHECK(xx.get() != NULL);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_stateL95_interpolation_constructor) {
    lorenz95::StateL95 xx1(*resol_, *vars_, *time_);
    boost::scoped_ptr<lorenz95::StateL95> xx2(new lorenz95::StateL95(*resol_, xx1));
    BOOST_CHECK(xx2.get() != NULL);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_stateL95_copy_constructor) {
    lorenz95::StateL95 xx(*resol_, *vars_, *time_);
    boost::scoped_ptr<lorenz95::StateL95> xx2(new lorenz95::StateL95(xx));
    BOOST_CHECK(xx2.get() != NULL);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_stateL95_destructor) {
    // nothing to test
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_stateL95_getField) {
    lorenz95::StateL95 xx(*resol_, *vars_, *time_);

    // there are 2 values in FieldL95: the 1st is the *resol_ value,
    // the 2nd is a vector of doubles initialised to 0.0, the size of the
    // vector is the *resol_ value
    BOOST_CHECK_EQUAL(xx.getField().resol(), resol_->npoints());
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_stateL95_validTime) {
    lorenz95::StateL95 xx(*resol_, *vars_, *time_);
    BOOST_CHECK_EQUAL(xx.validTime().toString(), date_str_);
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_stateL95_assignment) {
    util::DateTime tt(file_->getString("date"));
    lorenz95::StateL95 xx1(*resol_, *vars_, tt);
    xx1.read(*file_);

    // construct the second StateL95 object
    std::string date_string2("2014-09-14T09:35:00Z");
    util::DateTime dt2(date_string2);
    lorenz95::StateL95 xx2(*resol_, *vars_, dt2);

    xx2 = xx1;

    // initially, xx1 held data, xx2 was initialised to 0, now
    // ensure the two stateL95s are the same
    BOOST_CHECK_EQUAL(xx1.validTime(), xx2.validTime());
    BOOST_CHECK_EQUAL(xx1.getField().resol(), xx2.getField().resol());
    for (int i = 0; i < xx1.getField().resol(); ++i) {
      BOOST_CHECK_EQUAL(xx1.getField()[i], xx2.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  /*
  BOOST_AUTO_TEST_CASE(test_stateL95_interpolate) {
    lorenz95::StateL95 xx(*resol_, *vars_, *time_);

    int maxVecSize = 5;

    // construct the LocsL95 object
    std::vector<double> doubleVec(maxVecSize);
    // populate the vector with values from 0.1 to 0.5
    for(int i = 0; i < maxVecSize; ++i) {
      doubleVec[i] = ((i + 1) * 0.1);
    }
    lorenz95::LocsL95 locsL95(doubleVec);

    // construct the GomL95 object
    std::vector<int> intVec(maxVecSize);
    lorenz95::GomL95 gomL95(intVec);
    // populate the locval_ vector with values from 1.1 to 0.5
    for(int i = 0; i < intVec.size(); ++i) {
      gomL95[i] = ((i + 1) * 0.1);
    }

    int origCurrent = gomL95.current();

    xx.interpolate(locsL95, gomL95);

    for(int i = 0; i < locsL95->nobs(); ++i) {
      std::cout << "PMC: gomL95   " << gomL95[origCurrent + i] << std::endl;
      std::cout << "PMC: xx " << xx.getField()[i] << std::endl;
      BOOST_CHECK_EQUAL(gomL95[origCurrent + i], xx.getField()[i]);
    }
    BOOST_CHECK_EQUAL(gomL95.current(), (origCurrent + locsL95.nobs()));
    std::cout << "PMC: current()   " << gomL95.current() << std::endl;
    std::cout << "PMC: origCurrent " << origCurrent << std::endl;
    std::cout << "PMC: locs-nobs() " << locsL95.nobs() << std::endl;
  }
  */
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_stateL95_compound_assignment) {
    util::DateTime tt(file_->getString("date"));
    lorenz95::StateL95 xx(*resol_, *vars_, tt);
    xx.read(*file_);

    // construct the incrementL95 object
    lorenz95::IncrementL95 dx(*resol_, *vars_, tt);
    dx.read(*file_);

    xx += dx;

    // both xx and dx were initialised with the same data,
    // so when the two are added together xx should be double incL95
    for (int i = 0; i < xx.getField().resol(); ++i) {
      BOOST_CHECK_EQUAL(xx.getField()[i], 2.0 * dx.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_stateL95_read) {
    util::DateTime tt(file_->getString("date"));
    lorenz95::StateL95 xx(*resol_, *vars_, tt);
    xx.read(*file_);

    // to verify the information has been read correctly, we need to open
    // and read the file using ifstream functionality
    const std::string filename(file_->getString("filename"));
    std::ifstream inStream(filename.c_str());
    if (!inStream.is_open()) {
      BOOST_ERROR("read functionality cannot be determined");
    }

    int resolInt;
    inStream >> resolInt;

    std::string time;
    inStream >> time;

    std::vector<double> doubleVec(resolInt);
    for (int i = 0; i < resolInt; ++i) {
      inStream >> doubleVec[i];
    }
    inStream.close();

    for (int i = 0; i < resol_->npoints(); ++i) {
      BOOST_CHECK_EQUAL((xx.getField())[i], doubleVec[i]);
    }
  }
// -----------------------------------------------------------------------------
/*
  BOOST_AUTO_TEST_CASE(test_stateL95_stream_output) {
    lorenz95::StateL95 xx(*resol_, *vars_, *time_);
    xx.read(*file_);

    // use the operator<< method to write the value to a file

    std::filebuf fb;
    std::string filename("StateL95Test.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << xx;
    fb.close();

    std::vector<double> doubleVec(xx.getField().resol());
    for (int i = 0; i < xx.getField().resol(); ++i) {
      doubleVec[i] = xx.getField()[i];
    }
    std::vector<double>::iterator iter;
    iter = std::min_element(doubleVec.begin(), doubleVec.end());
    double min = *iter;
    iter = std::max_element(doubleVec.begin(), doubleVec.end());
    double max = *iter;
//    double avg = 0.0;

    // then read the values that were written to the file
    int idxStart;
    int idxEnd;
    double inputDouble = 0.0;
    std::string input;
    std::string inputTest;
    std::string inputTimeTest(" Valid time: " + file_->getString("date"));
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, input);  // ignore the first (blank) line
      getline(inputFile, input);

      BOOST_CHECK_EQUAL(input, inputTimeTest);

      getline(inputFile, input);
      // test min value
      idxStart = input.find("=");
      idxEnd = input.find(",");
      inputTest = input.substr(idxStart + 1, (idxEnd - 1) - idxStart);
      inputDouble = boost::lexical_cast<double>(inputTest);
      BOOST_CHECK_CLOSE(inputDouble, min, 0.0001);
      // test max value
      idxStart = input.find("=", idxEnd);
      idxEnd = input.find(",", idxEnd + 1);
      inputTest = input.substr(idxStart + 1, (idxEnd - 1) - idxStart);
      inputDouble = boost::lexical_cast<double>(inputTest);
      BOOST_CHECK_CLOSE(inputDouble, max, 0.0001);
      // test avg value
//      idxStart = input.find("=", idxEnd + 1);
//      idxEnd = input.find(",", idxEnd + 1);
//      inputTest = input.substr(idxStart + 1);
//      inputDouble = boost::lexical_cast<double>(inputTest);
//      BOOST_CHECK_CLOSE(inputDouble, avg, 0.0001);
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      BOOST_ERROR("operator<< functionality cannot be determined");
    }
    inputFile.close();
  }
*/
// -----------------------------------------------------------------------------
//  BOOST_AUTO_TEST_CASE(test_stateL95_step_traj) {
//    lorenz95::StateL95 xx(*resol_, *vars_, *time_);
//    xx.read(*file_);
//
//    // construct a ModelL95 object
//    eckit::LocalConfiguration modelCfg(TestConfig::config(), "model");
//    lorenz95::ModelL95 modelL95(*resol_, modelCfg);
//
//    // construct a ModelBias object
//    eckit::LocalConfiguration biasCfg(TestConfig::config(), "ModelBias");
//    lorenz95::ModelBias modelBias(*resol_, biasCfg);
//
//    // construct a ModelTrajectory object
//    lorenz95::ModelTrajectory modelTraj(true);
//    // copy construct a FieldL95 object to set params
//    // in the ModelTrajectory object
//    lorenz95::FieldL95 fieldL95(xx.getField());
//    modelTraj.set(fieldL95);
//
//    xx.stepTraj(modelL95, modelBias, modelTraj);
//
//    // create test data
//    modelL95.stepRK(fieldL95, modelBias, modelTraj);
//
//    // test xx data against newly created fieldL95 test data
//    for(int i = 0; i < xx.getField().resol(); ++i) {
//      BOOST_CHECK_EQUAL(xx.getField()[i], fieldL95[i]);
//    }
//  }
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
}  // namespace test
