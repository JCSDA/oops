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

  oops::instantiateObsErrorFactory<OBS>();

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
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

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    std::unique_ptr<Covar_> R(
      oops::ObsErrorFactory<OBS>::create(rconf, Test_::obspace()[jj]));

    // RMSE should be equal to the rms that was read from the file
    EXPECT(oops::is_close(R->getRMSE(), obserr.rms(), 1.e-10));

    // create random vector dy and its copies dy1, dy2
    ObsVector_ dy(Test_::obspace()[jj]);
    R->randomize(dy);
    ObsVector_ dy1(dy);
    ObsVector_ dy2(dy);
    oops::Log::test() << "Random vector dy: " << dy << std::endl;

    R->multiply(dy1);
    oops::Log::test() << "R*dy: " << dy1 << std::endl;
    R->inverseMultiply(dy1);
    // dy1 = R^{-1}*R*dy
    oops::Log::test() << "R^{-1}*R*dy: " << dy1 << std::endl;
    EXPECT(oops::is_close(dy1.rms(), dy.rms(), 1.e-10));

    R->inverseMultiply(dy2);
    oops::Log::test() << "R^{-1}*dy: " << dy2 << std::endl;
    R->multiply(dy2);
    // dy2 = R*R^P-1}*dy
    oops::Log::test() << "R*R^{-1}*dy: " << dy2 << std::endl;
    EXPECT(oops::is_close(dy2.rms(), dy.rms(), 1.e-10));
  }
}

// -----------------------------------------------------------------------------
/// Tests that the methods obserrors(), inverseVariance update() and save()
/// do what is expected.
template <typename OBS> void testAccessors() {
  typedef ObsTestsFixture<OBS>     Test_;
  typedef oops::ObsErrorBase<OBS>  Covar_;
  typedef oops::ObsVector<OBS>     ObsVector_;

  oops::instantiateObsErrorFactory<OBS>();

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ obserr(Test_::obspace()[jj], "ObsError");

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    std::unique_ptr<Covar_> R(
      oops::ObsErrorFactory<OBS>::create(rconf, Test_::obspace()[jj]));

    ObsVector_ dy(R->obserrors());
    oops::Log::test() << "ObsError: " << dy << std::endl;
    EXPECT(oops::is_close(dy.rms(), obserr.rms(), 1.e-10));

    ObsVector_ dy1(R->inverseVariance());
    oops::Log::test() << "inverseVariance: " << dy1 << std::endl;
    dy *= dy;
    dy.invert();
    EXPECT(oops::is_close(dy.rms(), dy1.rms(), 1.e-10));

    dy.ones();
    R->update(dy);
    oops::Log::test() << "R filled with ones: " << R->obserrors() << std::endl;
    EXPECT(oops::is_close(R->obserrors().rms(), R->inverseVariance().rms(), 1.e-10));
    EXPECT(oops::is_close(R->obserrors().rms(), dy.rms(), 1.e-10));

    R->save("Ones");
    ObsVector_ testOnes(Test_::obspace()[jj], "Ones");
    EXPECT(oops::is_close(dy.rms(), testOnes.rms(), 1.e-10));
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
    ts.emplace_back(CASE("interface/ObsErrorCovariance/testAccessors")
      { testAccessors<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_BASE_OBSERRORCOVARIANCE_H_
