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
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/FieldL95.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"
#include "oops/base/Variables.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
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
  std::unique_ptr<const eckit::LocalConfiguration> file_;
  std::unique_ptr<lorenz95::Resolution> resol_;
  std::string date_str_;
  std::unique_ptr<util::DateTime> time_;
  std::unique_ptr<oops::Variables> vars_;
};
// -----------------------------------------------------------------------------
CASE("test_StateL95") {
  StateTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_stateL95_constructor") {
    std::unique_ptr<lorenz95::StateL95>
      xx(new lorenz95::StateL95(*fix.resol_, *fix.vars_, *fix.time_));
    EXPECT(xx.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_stateL95_readin_constructor") {
    std::unique_ptr<lorenz95::StateL95>
      xx(new lorenz95::StateL95(*fix.resol_, *fix.vars_, *fix.file_));
    EXPECT(xx.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_stateL95_interpolation_constructor") {
    lorenz95::StateL95 xx1(*fix.resol_, *fix.vars_, *fix.time_);
    std::unique_ptr<lorenz95::StateL95> xx2(new lorenz95::StateL95(*fix.resol_, xx1));
    EXPECT(xx2.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_stateL95_copy_constructor") {
    lorenz95::StateL95 xx(*fix.resol_, *fix.vars_, *fix.time_);
    std::unique_ptr<lorenz95::StateL95> xx2(new lorenz95::StateL95(xx));
    EXPECT(xx2.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_stateL95_destructor") {
    // nothing to test
  }
// -----------------------------------------------------------------------------
  SECTION("test_stateL95_getField") {
    lorenz95::StateL95 xx(*fix.resol_, *fix.vars_, *fix.time_);

    // there are 2 values in FieldL95: the 1st is the *resol_ value,
    // the 2nd is a vector of doubles initialised to 0.0, the size of the
    // vector is the *resol_ value
    EXPECT(xx.getField().resol() == fix.resol_->npoints());
  }
// -----------------------------------------------------------------------------
  SECTION("test_stateL95_validTime") {
    lorenz95::StateL95 xx(*fix.resol_, *fix.vars_, *fix.time_);
    EXPECT(xx.validTime().toString() == fix.date_str_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_stateL95_assignment") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::StateL95 xx1(*fix.resol_, *fix.vars_, tt);
    xx1.read(*fix.file_);

    // construct the second StateL95 object
    std::string date_string2("2014-09-14T09:35:00Z");
    util::DateTime dt2(date_string2);
    lorenz95::StateL95 xx2(*fix.resol_, *fix.vars_, dt2);

    xx2 = xx1;

    // initially, xx1 held data, xx2 was initialised to 0, now
    // ensure the two stateL95s are the same
    EXPECT(xx1.validTime() == xx2.validTime());
    EXPECT(xx1.getField().resol() == xx2.getField().resol());
    for (int i = 0; i < xx1.getField().resol(); ++i) {
      EXPECT(xx1.getField()[i] == xx2.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_stateL95_compound_assignment") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::StateL95 xx(*fix.resol_, *fix.vars_, tt);
    xx.read(*fix.file_);

    // construct the incrementL95 object
    lorenz95::IncrementL95 dx(*fix.resol_, *fix.vars_, tt);
    dx.read(*fix.file_);

    xx += dx;

    // both xx and dx were initialised with the same data,
    // so when the two are added together xx should be double incL95
    for (int i = 0; i < xx.getField().resol(); ++i) {
      EXPECT(xx.getField()[i] == 2.0 * dx.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_stateL95_read") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::StateL95 xx(*fix.resol_, *fix.vars_, tt);
    xx.read(*fix.file_);

    // to verify the information has been read correctly, we need to open
    // and read the file using ifstream functionality
    const std::string filename(fix.file_->getString("filename"));
    std::ifstream inStream(filename.c_str());
    if (!inStream.is_open()) {
      oops::Log::error() << "read functionality cannot be determined" << std::endl;
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

    for (int i = 0; i < fix.resol_->npoints(); ++i) {
      EXPECT((xx.getField())[i] == doubleVec[i]);
    }
  }
// -----------------------------------------------------------------------------
/*
  SECTION("test_stateL95_stream_output") {
    lorenz95::StateL95 xx(*fix.resol_, *fix.vars_, *fix.time_);
    xx.read(*fix.file_);

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
    std::string inputTimeTest(" Valid time: " + fix.file_->getString("date"));
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, input);  // ignore the first (blank) line
      getline(inputFile, input);

      EXPECT(input == inputTimeTest);

      getline(inputFile, input);
      // test min value
      idxStart = input.find("=");
      idxEnd = input.find(",");
      inputTest = input.substr(idxStart + 1, (idxEnd - 1) - idxStart);
      inputDouble = boost::lexical_cast<double>(inputTest);
      EXPECT(oops::is_close(inputDouble, min, 0.000001));
      // test max value
      idxStart = input.find("=", idxEnd);
      idxEnd = input.find(",", idxEnd + 1);
      inputTest = input.substr(idxStart + 1, (idxEnd - 1) - idxStart);
      inputDouble = boost::lexical_cast<double>(inputTest);
      EXPECT(oops::is_close(inputDouble, max, 0.000001));
      // test avg value
//      idxStart = input.find("=", idxEnd + 1);
//      idxEnd = input.find(",", idxEnd + 1);
//      inputTest = input.substr(idxStart + 1);
//      inputDouble = boost::lexical_cast<double>(inputTest);
//      EXPECT(oops::is_close(inputDouble, avg, 0.000001));
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      oops::Log::error() << "operator<< functionality cannot be determined" << std::endl;
    }
    inputFile.close();
  }
*/
// -----------------------------------------------------------------------------
//  SECTION("test_stateL95_step_traj") {
//    lorenz95::StateL95 xx(*fix.resol_, *fix.vars_, *fix.time_);
//    xx.read(*fix.file_);
//
//    // construct a ModelL95 object
//    eckit::LocalConfiguration modelCfg(TestConfig::config(), "model");
//    lorenz95::ModelL95 modelL95(*fix.resol_, modelCfg);
//
//    // construct a ModelBias object
//    eckit::LocalConfiguration biasCfg(TestConfig::config(), "ModelBias");
//    lorenz95::ModelBias modelBias(*fix.resol_, biasCfg);
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
//      EXPECT(xx.getField()[i] == fieldL95[i]);
//    }
//  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
