/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_BASE_DEPARTURES_H_
#define TEST_BASE_DEPARTURES_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Departures.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename OBS> void testDepartures() {
  typedef ObsTestsFixture<OBS>     Test_;
  typedef oops::Departures<OBS>    Departures_;

  Departures_ y(Test_::obspace());
  Departures_ y2(Test_::obspace());
  Departures_ ydiff(Test_::obspace());

  // check assign and "-=" operator
  y.random();
  double ref_rms = y.rms();
  oops::Log::test() << "yrandom.rms=" << ref_rms << std::endl;

  y2.zero();
  oops::Log::test() << "yzero.rms=" << y2.rms() << std::endl;

  ydiff = y;
  ydiff -= y2;
  oops::Log::test() << "rms(yrandom-yzero)=" << ydiff.rms() << std::endl;

  EXPECT(ydiff.rms() != 0.0);

  // check pack operator
  Eigen::MatrixXd ypack = y.packEigen();
  oops::Log::info() << "ypack: " << ypack << std::endl;
  double rms = ypack.norm() / sqrt(ypack.size());
  oops::Log::test() << "rms(ypack)=" << rms << std::endl;
  oops::Log::test() << "y.rms()=" << ref_rms << std::endl;
  oops::Log::test() << "y.rms()-rms(ypack)" << ref_rms - rms << std::endl;

  EXPECT(oops::is_close(ref_rms, rms, 1.e-14));
}

template <typename OBS> class Departures : public oops::Test {
 public:
  Departures() {}
  virtual ~Departures() {}
 private:
  std::string testid() const {return "test::Departures<" + OBS::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();
    ts.emplace_back(CASE("base/Departures/testDepartures")
      { testDepartures<OBS>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_BASE_DEPARTURES_H_
