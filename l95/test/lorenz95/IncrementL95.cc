/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>  //  for std::unique_ptr

#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"
#include "oops/base/Variables.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class IncrementTestFixture : TestFixture {
 public:
  IncrementTestFixture() {
    file_.reset(new eckit::LocalConfiguration(TestConfig::config(), "state"));
    eckit::LocalConfiguration res(TestConfig::config(), "geometry");
    resol_.reset(new lorenz95::Resolution(res, oops::mpi::comm()));
    date_str_ = file_->getString("date");
    time_.reset(new util::DateTime(date_str_));
    vars_.reset(new oops::Variables());
  }
  ~IncrementTestFixture() {}
  std::unique_ptr<const eckit::LocalConfiguration> file_;
  std::unique_ptr<lorenz95::Resolution> resol_;
  std::string date_str_;
  std::unique_ptr<util::DateTime> time_;
  std::unique_ptr<oops::Variables> vars_;
};
// -----------------------------------------------------------------------------
CASE("test_IncrementL95") {
  IncrementTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_constructor") {
    std::unique_ptr<lorenz95::IncrementL95>
      dx(new lorenz95::IncrementL95(*fix.resol_, *fix.vars_, *fix.time_));
    EXPECT(dx.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_interpolation_constructor") {
    std::unique_ptr<lorenz95::IncrementL95>
      dx1(new lorenz95::IncrementL95(*fix.resol_, *fix.vars_, *fix.time_));
    std::unique_ptr<lorenz95::IncrementL95>
      dx2(new lorenz95::IncrementL95(*fix.resol_, *dx1));
    EXPECT(dx2.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_copy_constructor") {
    std::unique_ptr<lorenz95::IncrementL95>
      dx1(new lorenz95::IncrementL95(*fix.resol_, *fix.vars_, *fix.time_));
    std::unique_ptr<lorenz95::IncrementL95> dx2(new lorenz95::IncrementL95(*dx1));
    EXPECT(dx2.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_diff") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx(*fix.resol_, *fix.vars_, tt);
    dx.read(*fix.file_);

    // construct the first stateL95 object
    lorenz95::StateL95 xx1(*fix.resol_, *fix.vars_, tt);

    // read in the state config info
    xx1.read(*fix.file_);
    lorenz95::StateL95 xx2(xx1);

    // to vary the results a little, change the second StateL95 field values
    double fact = 0.75;
    dx *= fact;
    xx2 += dx;

    dx.diff(xx1, xx2);

    // to verify the diff method has worked correctly, we need
    // to open and read the file containing the FieldL95
    std::string filename(fix.file_->getString("filename"));
    std::ifstream inStream(filename.c_str());
    if (!inStream.is_open()) {
      oops::Log::error() << "diff functionality cannot be determined" << std::endl;
    }

    // we read in these two values but do not use them
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
      EXPECT(oops::is_close((dx.getField())[i],
                         doubleVec[i] - (doubleVec[i] + (doubleVec[i] * fact)),
                         1.0e-8));
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_zero") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx(*fix.resol_, *fix.vars_, tt);
    dx.read(*fix.file_);

    // first check that we have good data, ie, at least one element is non-zero
    bool goodData = false;
    for (int i = 0; i < fix.resol_->npoints() && goodData == false; ++i) {
      if ((dx.getField())[i] != 0) {
        goodData = true;
      }
    }

    if (!goodData) {
      oops::Log::error() <<
        "unable to test zero method, since test data is already all zero" << std::endl;;
    } else {
      dx.zero();

      for (int i = 0; i < fix.resol_->npoints(); ++i) {
        EXPECT(dx.getField()[i] == 0);
      }
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_zero_set_datetime") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx(*fix.resol_, *fix.vars_, tt);
    dx.read(*fix.file_);

    // first check that we have good data, ie, at least one element is non-zero
    bool goodData = false;
    for (int i = 0; i < fix.resol_->npoints() && goodData == false; ++i) {
      if (dx.getField()[i] != 0) {
        goodData = true;
      }
    }

    if (!goodData) {
      oops::Log::error() <<
        "unable to test zero method, since test data is already all zero" << std::endl;
    } else {
      const std::string modified_date_string("2010-01-01T10:35:00Z");
      const util::DateTime dtModified(modified_date_string);
      dx.zero(dtModified);

      for (int i = 0; i < fix.resol_->npoints(); ++i) {
        EXPECT(dx.getField()[i] == 0);
      }

      EXPECT(dx.validTime().toString() != fix.date_str_);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_assignment") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx1(*fix.resol_, *fix.vars_, tt);
    dx1.read(*fix.file_);

    // construct the second dx object
    lorenz95::IncrementL95 dx2(*fix.resol_, *fix.vars_, tt);
    double fact = 0.75;
    dx2.read(*fix.file_);
    dx2 *= fact;

    dx1 = dx2;

    for (int i = 0; i < dx1.getField().resol(); ++i) {
      EXPECT(dx1.getField()[i] == dx2.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_compound_assignment_add") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx1(*fix.resol_, *fix.vars_, tt);
    dx1.read(*fix.file_);

    // copy construct the second stateL95 object
    lorenz95::IncrementL95 dx2(dx1);

    dx1 += dx2;

    // since the two IncrementL95 objects started off with the same data,
    // once they've been added together incL591 will be double what incL952 is
    for (int i = 0; i < dx1.getField().resol(); ++i) {
      EXPECT(dx1.getField()[i] == 2.0 * dx2.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_compound_assignment_subtract") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx1(*fix.resol_, *fix.vars_, tt);
    dx1.read(*fix.file_);

    // copy construct the second stateL95 object
    lorenz95::IncrementL95 dx2(dx1);

    dx1 -= dx2;

    // since the two IncrementL95 objects started off with the same data,
    // once incL952 has been subtracted from incL951, the result is zero
    for (int i = 0; i < dx1.getField().resol(); ++i) {
      EXPECT(dx1.getField()[i] == 0.0);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_compound_assignment_multiply") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx(*fix.resol_, *fix.vars_, tt);
    dx.read(*fix.file_);

    // create a copy of the original data for testing against
    std::vector<double> testData(dx.getField().resol());
    for (unsigned int ii = 0; ii < testData.size(); ++ii) {
      testData.at(ii) = dx.getField()[ii];
    }

    double fact = 0.75;
    dx *= fact;

    for (int ii = 0; ii < dx.getField().resol(); ++ii) {
      EXPECT(dx.getField()[ii] == testData.at(ii) * fact);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_axpy") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx1(*fix.resol_, *fix.vars_, tt);
    dx1.read(*fix.file_);

    // copy construct the second stateL95 object
    lorenz95::IncrementL95 dx2(dx1);

    double fact = 0.75;
    dx1.axpy(fact, dx2);

    for (int i = 0; i < dx1.getField().resol(); ++i) {
      EXPECT(dx1.getField()[i] == dx2.getField()[i] + fact * dx2.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_dot_product_with") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx1(*fix.resol_, *fix.vars_, tt);
    dx1.read(*fix.file_);

    // copy construct the second stateL95 object
    lorenz95::IncrementL95 dx2(dx1);

    double dpwResult = dx1.dot_product_with(dx2);

    // prepare a value to test against
    double testResult = 0.0;
    for (int i = 0; i < dx1.getField().resol(); ++i) {
      testResult += (dx1.getField()[i] * dx2.getField()[i]);
    }

    EXPECT(dpwResult == testResult);
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_schur_product_with") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx1(*fix.resol_, *fix.vars_, tt);
    dx1.read(*fix.file_);

    // copy construct the second stateL95 object
    lorenz95::IncrementL95 dx2(dx1);

    dx1.schur_product_with(dx2);

    // both incL951 and incL952 started off with the same data in x_,
    // so to test incL951 against incL952xincL952 is a valid test
    for (int i = 0; i < dx1.getField().resol(); ++i) {
      EXPECT(dx1.getField()[i] == dx2.getField()[i] * dx2.getField()[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_read") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx(*fix.resol_, *fix.vars_, tt);
    dx.read(*fix.file_);

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
      EXPECT(dx.getField()[i] == doubleVec[i]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_write") {
    util::DateTime tt(fix.file_->getString("date"));
    lorenz95::IncrementL95 dx(*fix.resol_, *fix.vars_, tt);
    dx.read(*fix.file_);

    eckit::LocalConfiguration opFileCfg(TestConfig::config(), "output file");
    dx.write(opFileCfg);

    // Should read back in and compare values
  }
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_validTime") {
    lorenz95::IncrementL95 dx(*fix.resol_, *fix.vars_, *fix.time_);
    EXPECT(dx.validTime().toString() == fix.date_str_);
  }
// -----------------------------------------------------------------------------
/*
  SECTION("test_incrementL95_stream_output") {
    lorenz95::IncrementL95 dx(*fix.resol_, *fix.vars_, *fix.time_);

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("IncrementL95Test.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << dx;
    fb.close();

    // then read the value that was written to the file
    std::string input;
    std::string inputTest(" Valid time: " + fix.date_str_);
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, input);  // ignore the first (blank) line
      getline(inputFile, input);

      EXPECT(input == inputTest);
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      oops::Log::error() << "operator<< functionality cannot be determined" << std::endl;
    }
    inputFile.close();
  }
*/
// -----------------------------------------------------------------------------
  SECTION("test_incrementL95_getField") {
    lorenz95::IncrementL95 dx(*fix.resol_, *fix.vars_, *fix.time_);

    // there are 2 values in FieldL95: the 1st is the *resol_ value,
    // the 2nd is a vector of doubles initialised to 0.0, the size of the
    // vector is the *resol_ value (we're just checking the final one)
    EXPECT(dx.getField().resol() == fix.resol_->npoints());
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
