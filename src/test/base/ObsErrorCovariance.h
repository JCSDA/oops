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

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/exception/Exceptions.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/ObsError.h"
#include "oops/runs/Test.h"
#include "oops/util/missingValues.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// Tests creation and destruction of ObsErrorCovariances
template <typename OBS> void testConstructor() {
  typedef ObsTestsFixture<OBS>                 Test_;
  typedef oops::ObsError<OBS>                  Covar_;

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    std::unique_ptr<Covar_> R = std::make_unique<Covar_>(rconf,
                                                         Test_::obspace()[jj]);
    EXPECT(R.get());
    oops::Log::info() << "Testing ObsError: " << *R << std::endl;
    R.reset();
    EXPECT(!R.get());
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/// Test whether correlation matrix is being read correctly:
/// multiply R by unit obs vector and compare to reference vector
template <typename OBS> void testReader() {
  typedef ObsTestsFixture<OBS>                 Test_;
  typedef oops::ObsError<OBS>                  Covar_;
  typedef oops::ObsVector<OBS>                 ObsVector_;

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    if (!conf[jj].has("obs error test")) {
        const std::string name = Test_::obspace()[jj].obsname();
        oops::Log::info() << name + ": Test Reader not found" << std::endl;
        continue;
    }
    const eckit::LocalConfiguration testConf(conf[jj], "obs error test");

      if (testConf.getBool("test reader", false)) {
        const eckit::LocalConfiguration rconf(conf[jj], "obs error");
        Covar_ R(rconf, Test_::obspace()[jj]);

        // Read in obs vector from Obsvalues , this will be a unit vector e.g [1, 0, 0]
        ObsVector_ unit(Test_::obspace()[jj], "ObsValue");
        R.multiply(unit);
        ObsVector_ mask(Test_::obspace()[jj]);
        mask.zero();
        std::vector<double> unitVec;
        unit.maskAndSerialize(mask, unitVec);
        oops::Log::info() << "Column of R matrix: " << std::endl << unitVec << std::endl;
        std::vector<float> refVec = testConf.getFloatVector("reference");

        // unitVec after Multiplication with R should equal reference vector
        for (size_t i = 0; i < unitVec.size(); i++) {
            EXPECT_EQUAL(unitVec[i], refVec[i]);
        }
    }
  }
}

// -----------------------------------------------------------------------------
/// Tests that \f$R*R^{-1}*dy = dy\f$ and \f$R^{-1}*R*dy = dy\f$
template <typename OBS> void testMultiplies() {
  typedef ObsTestsFixture<OBS>                 Test_;
  typedef oops::ObsError<OBS>                  Covar_;
  typedef oops::ObsVector<OBS>                 ObsVector_;

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ obserr(Test_::obspace()[jj], "ObsError");

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    Covar_ R(rconf, Test_::obspace()[jj]);

    double RMSE_tolerance = rconf.getDouble("Obs Error test tolerance", 1.e-10);
    // RMSE should be equal to the rms that was read from the file
    EXPECT(oops::is_close(R.getRMSE(), obserr.rms(), RMSE_tolerance));

    // create random vector dy and its copies dy1, dy2
    ObsVector_ dy(Test_::obspace()[jj], "");
    R.randomize(dy);
    ObsVector_ dy1(dy);
    ObsVector_ dy2(dy);
    oops::Log::info() << "Random vector dy: " << dy << std::endl;

    R.multiply(dy1);
    oops::Log::info() << "R*dy: " << dy1 << std::endl;
    R.inverseMultiply(dy1);
    // dy1 = R^{-1}*R*dy
    oops::Log::info() << "R^{-1}*R*dy: " << dy1 << std::endl;
    EXPECT(oops::is_close(dy1.rms(), dy.rms(), RMSE_tolerance));

    R.inverseMultiply(dy2);
    oops::Log::info() << "R^{-1}*dy: " << dy2 << std::endl;
    R.multiply(dy2);
    // dy2 = R*R^{-1}*dy
    oops::Log::info() << "R*R^{-1}*dy: " << dy2 << std::endl;
    EXPECT(oops::is_close(dy2.rms(), dy.rms(), RMSE_tolerance));
  }
}

// -----------------------------------------------------------------------------
/// Tests that the methods obserrors(), inverseVariance update() and save()
/// do what is expected.
template <typename OBS> void testAccessors() {
  typedef ObsTestsFixture<OBS>                 Test_;
  typedef oops::ObsError<OBS>                  Covar_;
  typedef oops::ObsVector<OBS>                 ObsVector_;

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ obserr(Test_::obspace()[jj], "ObsError");

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    double RMSE_tolerance = rconf.getDouble("RMSE tolerance", 1.e-10);
    Covar_ R(rconf, Test_::obspace()[jj]);

    ObsVector_ dy(R.obserrors());
    oops::Log::info() << "ObsError: " << dy << std::endl;
    EXPECT(oops::is_close(dy.rms(), obserr.rms(), RMSE_tolerance));

    ObsVector_ dy1(R.inverseVariance());
    oops::Log::info() << "inverseVariance: " << dy1 << std::endl;
    dy *= dy;
    dy.invert();
    EXPECT(oops::is_close(dy.rms(), dy1.rms(), RMSE_tolerance));

    dy.ones();
    R.update(dy);
    oops::Log::info() << "R filled with ones: " << R.obserrors() << std::endl;
    EXPECT(oops::is_close(R.obserrors().rms(), R.inverseVariance().rms(), RMSE_tolerance));
    EXPECT(oops::is_close(R.obserrors().rms(), dy.rms(), RMSE_tolerance));

    R.save("Ones");
    ObsVector_ testOnes(Test_::obspace()[jj], "Ones");
    EXPECT(oops::is_close(dy.rms(), testOnes.rms(), RMSE_tolerance));
  }
}

// -----------------------------------------------------------------------------
/// Test localization of R.
template <typename OBS> void testLocalize() {
  typedef ObsTestsFixture<OBS>                 Test_;
  typedef oops::ObsError<OBS>                  Covar_;
  typedef oops::ObsVector<OBS>                 ObsVector_;

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  bool testLocalize;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    if (conf[jj].has("obs error test")) {
      const eckit::LocalConfiguration testConf(conf[jj], "obs error test");
      testLocalize = testConf.getBool("test localization", true);
    } else {
      testLocalize = true;
    }
    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    Covar_ R(rconf, Test_::obspace()[jj]);

    ObsVector_ maskmissing(Test_::obspace()[jj]);
    maskmissing.ones();
    maskmissing *= util::missingValue<double>();
    if (testLocalize) {
      R.localize(maskmissing);
      EXPECT(R.localDim() == 0);
    } else {
      EXPECT_THROWS(R.localize(maskmissing));
    }

    ObsVector_ maskones(Test_::obspace()[jj]);
    maskones.ones();
    if (testLocalize) {
      R.localize(maskones);
    } else {
      EXPECT_THROWS(R.localize(maskones));
    }

    ObsVector_ masknegative(Test_::obspace()[jj]);
    masknegative.ones();
    masknegative *= -1;
    if (testLocalize) {
        EXPECT_THROWS_AS(R.localize(masknegative), eckit::BadValue);
    } else {
        EXPECT_THROWS(R.localize(masknegative));
    }
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsErrorCovariance : public oops::Test {
  typedef ObsTestsFixture<OBS>     Test_;

 public:
  using oops::Test::Test;
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
    ts.emplace_back(CASE("interface/ObsErrorCovariance/testReader")
      { testReader<OBS>(); });
    ts.emplace_back(CASE("interface/ObsErrorCovariance/testLocalize")
      { testLocalize<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_BASE_OBSERRORCOVARIANCE_H_
