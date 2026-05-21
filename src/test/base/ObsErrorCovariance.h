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
#include <numeric>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/exception/Exceptions.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/ObsError.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
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

    // update/finalize R from dy_update (((after qc)))
    R.update(obserr);

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
        EXPECT_THROWS_AS_ANY_RANK(R.localize(masknegative), eckit::BadValue);
    } else {
        EXPECT_THROWS(R.localize(masknegative));
    }
  }
}

// -----------------------------------------------------------------------------
/// Test localInverseMultiply routine.
template <typename OBS> void testLocalInverseMultiply() {
  typedef ObsTestsFixture<OBS>                 Test_;
  typedef oops::ObsError<OBS>                  Covar_;
  typedef oops::ObsVector<OBS>                 ObsVector_;

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("observations", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    if (!conf[jj].has("obs error test")) {
        const std::string name = Test_::obspace()[jj].obsname();
        oops::Log::info() << name + ": obs error test not found" << std::endl;
        continue;
    }
    const eckit::LocalConfiguration testConf(conf[jj], "obs error test");

    // It is only possible to run this test if the `localize()` function has been
    // implemented for the obs error class.
    if (!testConf.getBool("test localization", false)) {
      continue;
    }

    const eckit::LocalConfiguration rconf(conf[jj], "obs error");
    Covar_ R(rconf, Test_::obspace()[jj]);

    // Perform localization.
    ObsVector_ maskones(Test_::obspace()[jj]);
    maskones.ones();
    R.localize(maskones);

    const bool testLocalInverseMultiply = testConf.has("RMS localInverseMultiply");
    // If testing localInverseMultiply has not been requested, expect that calling
    // `R.localInverseMultiply()` for an obs error class will throw an exception.
    // This ensures that the user will not forget to add a reference RMS
    // if they fill in the `localInverseMultiply()` function for the obs error class.
    if (!testLocalInverseMultiply) {
      EXPECT_THROWS(R.localInverseMultiply(Eigen::MatrixXf(1, R.localDim())));
      continue;
    }

    const float RMSref = testConf.getFloat("RMS localInverseMultiply");
    // A negative reference value of RMS indicates that the test should not be performed.
    // This is useful when dealing with artificial data sets that do not contain the
    // ObsValue and GsiHofX groups.
    if (RMSref < 0.0) {
      continue;
    }

    // Compute O-B.
    ObsVector_ omb(Test_::obspace()[jj], "ObsValue");
    const ObsVector_ bkg(Test_::obspace()[jj], "GsiHofX");
    omb -= bkg;

    // Obs Error used for masking.
    const ObsVector_ err(Test_::obspace()[jj], "ObsError");

    // Convert omb into its localized version (local_omb).
    std::vector<double> vals;
    omb.maskAndSerialize(err, vals);
    const Eigen::VectorXd local_omb = Eigen::Map<Eigen::VectorXd>(vals.data(), vals.size());

    // Convert local_omb into an Eigen::MatrixXf.
    const Eigen::MatrixXf local_omb_f = static_cast<Eigen::MatrixXf>(local_omb.cast<float>());

    // Run localInverseMultiply.
    const Eigen::MatrixXf LIM = R.localInverseMultiply(local_omb_f.transpose());

    // Expect the result of R^-1 OmB will have a single row and at least one column.
    EXPECT(LIM.rows() == 1);
    EXPECT(LIM.cols() > 0);

    // Compute RMS of R^-1 OmB.
    const Eigen::VectorXf LIM0 = static_cast<Eigen::VectorXf>(LIM.row(0));
    std::vector<float> LIMvector(LIM0.data(), LIM0.data() + LIM0.rows());
    oops::mpi::allGatherv(Test_::comm(), LIMvector);
    EXPECT(LIMvector.size() == err.nobs());
    const float RMS = std::sqrt(std::inner_product(LIMvector.begin(),
                                                   LIMvector.end(),
                                                   LIMvector.begin(),
                                                   0.0) / static_cast<float>(LIMvector.size()));

    // Compare RMS against reference value.
    const float RMStol = rconf.getFloat("Obs Error test tolerance", 1e-5);
    EXPECT(oops::is_close(RMS, RMSref, RMStol));
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
    ts.emplace_back(CASE("interface/ObsErrorCovariance/testLocalInverseMultiply")
      { testLocalInverseMultiply<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_BASE_OBSERRORCOVARIANCE_H_
