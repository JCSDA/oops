/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_OBSSPACE_H_
#define TEST_INTERFACE_OBSSPACE_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/ObsSpace.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// \brief tests constructor, window accessor methods and prints
template <typename OBS> void testConstructor() {
  typedef ObsTestsFixture<OBS> Test_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    oops::Log::info() << "Testing ObsSpace: " << Test_::obspace()[jj] << std::endl;
    EXPECT(Test_::obspace()[jj].windowStart() == Test_::timeWindow().start());
    EXPECT(Test_::obspace()[jj].windowEnd() == Test_::timeWindow().end());
  }
}

// -----------------------------------------------------------------------------
/// \brief tests that ObsSpaces created on subwindows have the same number obs as
///        ObsSpaces created on the whole window
template <typename OBS> void testSubwindows() {
  typedef ObsTestsFixture<OBS>   Test_;
  typedef oops::ObsSpace<OBS>    ObsSpace_;
  typedef oops::ObsVector<OBS>   ObsVector_;

  const util::DateTime tbgn = Test_::timeWindow().start();
  const util::DateTime tmid = Test_::timeWindow().midpoint();
  const util::DateTime tend = Test_::timeWindow().end();
  const util::TimeWindow subWin1 = Test_::timeWindow().createSubWindow(tbgn, tmid);
  const util::TimeWindow subWin2 = Test_::timeWindow().createSubWindow(tmid, tend);
  oops::Log::info() << "Testing subwindows: " << std::endl;
  oops::Log::info() << subWin1 << std::endl;
  oops::Log::info() << subWin2 << std::endl;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration obsconfig(Test_::config(jj), "obs space");
    ObsSpace_ obspace1(obsconfig, oops::mpi::world(), subWin1);
    ObsSpace_ obspace2(obsconfig, oops::mpi::world(), subWin2);

    /// Create ObsVectors for each of the ObsSpaces, to compare nobs
    ObsVector_ ovec(Test_::obspace()[jj]);
    ObsVector_ ovec1(obspace1);
    ObsVector_ ovec2(obspace2);
    oops::Log::info() << Test_::obspace()[jj].obsname() << " nobs(all): " << ovec.nobs()
                      << " nobs(1st subwindow): " << ovec1.nobs()
                      << " nobs(2nd subwindow): " << ovec2.nobs() << std::endl;
    EXPECT_EQUAL(ovec1.nobs() + ovec2.nobs(), ovec.nobs());
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> class ObsSpace : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;
 public:
  using oops::Test::Test;
  virtual ~ObsSpace() {}
 private:
  std::string testid() const override {return "test::ObsSpace<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsSpace/testConstructor")
      { testConstructor<OBS>(); });
    ts.emplace_back(CASE("interface/ObsSpace/testSubwindows")
      { testSubwindows<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_OBSSPACE_H_
