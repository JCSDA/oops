/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_BASE_OBSERRORCOVARIANCE_H_
#define TEST_BASE_OBSERRORCOVARIANCE_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/testing/Test.h"
#include "oops/base/ObsErrorBase.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/ObsVector.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// Tests creation and destruction of ObsErrorCovariances
template <typename OBS> void testConstructor() {
  typedef ObsTestsFixture<OBS>     Test_;
  typedef oops::ObsErrorBase<OBS>  Covar_;
  typedef oops::ObsVector<OBS>     ObsVector_;

  oops::instantiateObsErrorFactory<OBS>();

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ obserr(Test_::obspace()[jj], "ObsError");
    obserr.save("EffectiveError");

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    std::unique_ptr<Covar_> R(
      oops::ObsErrorFactory<OBS>::create(rconf, Test_::obspace()[jj]));
    EXPECT(R.get());
    oops::Log::test() << "Testing ObsError: " << *R << std::endl;
    R.reset();
    EXPECT(!R.get());
  }
}

// -----------------------------------------------------------------------------
/// Tests that \f$R*R^{-1}*dy = dy\f$ and \f$R^{-1}*R*dy = dy\f$
template <typename OBS> void testMultiplies() {
  typedef ObsTestsFixture<OBS>     Test_;
  typedef oops::ObsErrorBase<OBS>  Covar_;
  typedef oops::ObsVector<OBS>     ObsVector_;

  oops::instantiateObsErrorFactory<OBS>();

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ obserr(Test_::obspace()[jj], "ObsError");
    obserr.save("EffectiveError");

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    std::unique_ptr<Covar_> R(
      oops::ObsErrorFactory<OBS>::create(rconf, Test_::obspace()[jj]));

    // RMSE should be equal to the rms that was read from the file
    EXPECT(oops::is_close(R->getRMSE(), obserr.rms(), 1.e-10));

    // create random vector dy and its copies dy1, dy2
    ObsVector_ dy(Test_::obspace()[jj]);
    dy.random();
    ObsVector_ dy1(dy);
    ObsVector_ dy2(dy);
    oops::Log::info() << "Random vector dy: " << dy << std::endl;

    R->multiply(dy1);
    R->inverseMultiply(dy1);
    // dy1 = R^{-1}*R*dy
    oops::Log::info() << "R^{-1}*R*dy: " << dy1 << std::endl;
    EXPECT(oops::is_close(dy1.rms(), dy.rms(), 1.e-10));

    R->inverseMultiply(dy2);
    R->multiply(dy2);
    // dy2 = R*R^P-1}*dy
    oops::Log::info() << "R*R^{-1}*dy: " << dy2 << std::endl;
    EXPECT(oops::is_close(dy2.rms(), dy.rms(), 1.e-10));
  }
}


// -----------------------------------------------------------------------------

template <typename OBS>
class ObsErrorCovariance : public oops::Test {
  typedef ObsTestsFixture<OBS>     Test_;
 public:
  ObsErrorCovariance() {}
  virtual ~ObsErrorCovariance() {}
 private:
  std::string testid() const override {return "test::ObsErrorCovariance<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsErrorCovariance/testConstructor")
      { testConstructor<OBS>(); });
    ts.emplace_back(CASE("interface/ObsErrorCovariance/testMultiplies")
      { testMultiplies<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_BASE_OBSERRORCOVARIANCE_H_
