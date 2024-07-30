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

#include "eckit/testing/Test.h"
#include "oops/base/ObsError.h"
#include "oops/base/ObsVector.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

class ObsErrorTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorTestParameters, Parameters)
 public:
  oops::Parameter<bool> testReader{"test reader", false, this};
  oops::OptionalParameter<std::vector<float>> refVec{"reference", this};
};

// -----------------------------------------------------------------------------
/// Tests creation and destruction of ObsErrorCovariances
template <typename OBS> void testConstructor() {
  typedef ObsTestsFixture<OBS>                 Test_;
  typedef oops::ObsError<OBS>                  Covar_;
  typedef oops::ObsErrorParametersWrapper<OBS> Parameters_;

  oops::instantiateObsErrorFactory<OBS>();

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    Parameters_ rparams;
    rparams.validateAndDeserialize(rconf);
    std::unique_ptr<Covar_> R = std::make_unique<Covar_>(rparams.obsErrorParameters,
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
  typedef oops::ObsErrorParametersWrapper<OBS> Parameters_;
  typedef oops::ObsVector<OBS>                 ObsVector_;

  oops::instantiateObsErrorFactory<OBS>();

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    if (!conf[jj].has("obs error test")) {
        const std::string name = Test_::obspace()[jj].obsname();
        oops::Log::info() << name + ": Test Reader not found" << std::endl;
        continue;
    }
    const eckit::LocalConfiguration errconf(conf[jj], "obs error test");
    ObsErrorTestParameters testParams;
    testParams.validateAndDeserialize(errconf);

    if (testParams.testReader.value()) {
        const eckit::LocalConfiguration rconf(conf[jj], "obs error");
        Parameters_ rparams;
        rparams.validateAndDeserialize(rconf);
        Covar_ R(rparams.obsErrorParameters, Test_::obspace()[jj]);

        // Read in obs vector from Obsvalues , this will be a unit vector e.g [1, 0, 0]
        ObsVector_ unit(Test_::obspace()[jj], "ObsValue");
        R.multiply(unit);
        ObsVector_ mask(Test_::obspace()[jj]);
        mask.zero();
        Eigen::VectorXd unitVec = unit.packEigen(mask);
        oops::Log::info() << "Column of R matrix: " << std::endl << unitVec << std::endl;
        std::vector<float> refVec = testParams.refVec.value().value();

        // unitVec after Multiplication with R should equal reference vector
        for (int i = 0; i < unitVec.size(); i++) {
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
  typedef oops::ObsErrorParametersWrapper<OBS> Parameters_;
  typedef oops::ObsVector<OBS>                 ObsVector_;

  oops::instantiateObsErrorFactory<OBS>();

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ obserr(Test_::obspace()[jj], "ObsError");

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    Parameters_ rparams;
    rparams.validateAndDeserialize(rconf);
    Covar_ R(rparams.obsErrorParameters, Test_::obspace()[jj]);

    // RMSE should be equal to the rms that was read from the file
    EXPECT(oops::is_close(R.getRMSE(), obserr.rms(), 1.e-10));

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
    EXPECT(oops::is_close(dy1.rms(), dy.rms(), 1.e-10));

    R.inverseMultiply(dy2);
    oops::Log::info() << "R^{-1}*dy: " << dy2 << std::endl;
    R.multiply(dy2);
    // dy2 = R*R^P-1}*dy
    oops::Log::info() << "R*R^{-1}*dy: " << dy2 << std::endl;
    EXPECT(oops::is_close(dy2.rms(), dy.rms(), 1.e-10));
  }
}

// -----------------------------------------------------------------------------
/// Tests that the methods obserrors(), inverseVariance update() and save()
/// do what is expected.
template <typename OBS> void testAccessors() {
  typedef ObsTestsFixture<OBS>                 Test_;
  typedef oops::ObsError<OBS>                  Covar_;
  typedef oops::ObsErrorParametersWrapper<OBS> Parameters_;
  typedef oops::ObsVector<OBS>                 ObsVector_;

  oops::instantiateObsErrorFactory<OBS>();

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ obserr(Test_::obspace()[jj], "ObsError");

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    Parameters_ rparams;
    rparams.validateAndDeserialize(rconf);
    Covar_ R(rparams.obsErrorParameters, Test_::obspace()[jj]);

    ObsVector_ dy(R.obserrors());
    oops::Log::info() << "ObsError: " << dy << std::endl;
    EXPECT(oops::is_close(dy.rms(), obserr.rms(), 1.e-10));

    ObsVector_ dy1(R.inverseVariance());
    oops::Log::info() << "inverseVariance: " << dy1 << std::endl;
    dy *= dy;
    dy.invert();
    EXPECT(oops::is_close(dy.rms(), dy1.rms(), 1.e-10));

    dy.ones();
    R.update(dy);
    oops::Log::info() << "R filled with ones: " << R.obserrors() << std::endl;
    EXPECT(oops::is_close(R.obserrors().rms(), R.inverseVariance().rms(), 1.e-10));
    EXPECT(oops::is_close(R.obserrors().rms(), dy.rms(), 1.e-10));

    R.save("Ones");
    ObsVector_ testOnes(Test_::obspace()[jj], "Ones");
    EXPECT(oops::is_close(dy.rms(), testOnes.rms(), 1.e-10));
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
  }

  void clear() const override {
    Test_::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_BASE_OBSERRORCOVARIANCE_H_
