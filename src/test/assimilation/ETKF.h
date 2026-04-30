/*
 * (C) Crown copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_ASSIMILATION_ETKF_H_
#define TEST_ASSIMILATION_ETKF_H_

#include <Eigen/Dense>
#include <chrono>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/assimilation/ETKFLinearAlgebra.h"
#include "oops/assimilation/gletkfInterface.h"
#include "oops/assimilation/LETKFSolver.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"

namespace test {

  auto timeDifference = [](auto tA, auto tB) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(tB - tA).count();
  };

  class ObsErrorBogus {
    // Acts like a diagonal obs error matrix with variances sequentially increasing from 1 up to
    // nobs. Only public facing function is localInverseMultiply as that is all the linear
    // algebra tests use.
   public:
    explicit ObsErrorBogus(const int& nobs)
      : invVarR(Eigen::VectorXd::LinSpaced(nobs, 1.0, nobs)) {}
    Eigen::MatrixXd localInverseMultiply(const Eigen::MatrixXf & zz) const {
      Eigen::MatrixXf zzRinv(zz.rows(), zz.cols());
      for (int ii = 0; ii < zz.rows(); ++ii) {
        zzRinv(ii, Eigen::placeholders::all) = zz(ii, Eigen::placeholders::all)
        .cwiseProduct(invVarR.cast<float>().transpose());
      }
      return zzRinv.cast<double>();
    }
   private:
    const Eigen::VectorXd invVarR;
  };

  void compareWeights(const Eigen::VectorXf & wa, const Eigen::VectorXf & wa_f,
                      const Eigen::MatrixXf & Wa, const Eigen::MatrixXf & Wa_f) {
    const float wa_tol = 1.0e-6;
    const float Wa_tol = 1.0e-5;
    for (int i = 0; i < wa.rows(); ++i) {
      EXPECT(oops::is_close_absolute(wa(i), wa_f(i), wa_tol));
    }
    for (int i = 0; i < Wa.rows(); ++i) {
      for (int j = 0; j < Wa.cols(); ++j) {
        EXPECT(oops::is_close_absolute(Wa(i, j), Wa_f(i, j), Wa_tol));
      }
    }
  }

  void LETKF(const int nobs, const int nens) {
    srand(1);
    const bool singularValueDecomposition = false;

    const Eigen::VectorXd dy = Eigen::VectorXd::Random(nobs);
    const Eigen::MatrixXf Yb = Eigen::MatrixXf::Random(nens, nobs);
    const Eigen::VectorXd invVarR = Eigen::VectorXd::LinSpaced(nobs, 1.0, nobs);
    const ObsErrorBogus R(nobs);
    Eigen::VectorXf wa(nens);
    Eigen::MatrixXf Wa(nens, nens);

    // Eigen implementation

    const auto times = oops::detLETKF_computeWeights(
        dy.cast<float>(), Yb, R, (nens - 1) / 1.0,
        singularValueDecomposition, wa, Wa);

    // LAPACK implementation

    constexpr int neig = 1;
    constexpr int getkf_inflation = 0;
    constexpr int denkf = 0;
    constexpr int getkf = 0;
    constexpr float infl = 1.0;

    Eigen::VectorXf wa_f(nens);
    Eigen::MatrixXf Wa_f(nens, nens);

    const auto tL0 = std::chrono::system_clock::now();

    // cast eigen<double> to eigen<float>
    const Eigen::VectorXf dy_f = dy.cast<float>();
    const Eigen::VectorXf invVarR_f = invVarR.cast<float>();

    oops::letkf_core_f90(nobs,
                         Yb.data(),
                         Yb.data(),
                         dy_f.data(),
                         wa_f.data(),
                         Wa_f.data(),
                         invVarR_f.data(),
                         nens,
                         neig,
                         getkf_inflation,
                         denkf,
                         getkf,
                         infl);

    const auto tL1 = std::chrono::system_clock::now();

    // Timing information

    oops::Log::info() << std::endl;
    oops::Log::info() << "LETKF timing (ms)" << std::endl;
    oops::Log::info() << std::endl;
    oops::Log::info() << "Eigen:" << std::endl;
    oops::Log::info() << " compute Y^T R^-1 Y: " << timeDifference(times[0], times[1]) << std::endl;
    oops::Log::info() << " eigendecomposition: " << timeDifference(times[1], times[2]) << std::endl;
    oops::Log::info() << " compute Pa: " << timeDifference(times[2], times[3]) << std::endl;
    oops::Log::info() << " compute weights: " << timeDifference(times[3], times[4]) << std::endl;
    oops::Log::info() << "Eigen total: " << timeDifference(times[0], times[4]) << std::endl;
    oops::Log::info() << std::endl;
    oops::Log::info() << "LAPACK total: " << timeDifference(tL0, tL1) << std::endl;
    oops::Log::info() << std::endl;

    // Compare the weights produced by the two implementations

    compareWeights(wa, wa_f, Wa, Wa_f);
  }

  void GETKF(const int nobs, const int nens, const int neig) {
    srand(1);
    const bool svd = false;

    const int nana = nens * neig;
    constexpr float infl = 1.0;
    const Eigen::VectorXd dy = Eigen::VectorXd::Random(nobs);
    const Eigen::MatrixXf YbOrig = Eigen::MatrixXf::Random(nens, nobs);
    // modulated ensemble
    const Eigen::MatrixXf Yb = Eigen::MatrixXf::Random(nana, nobs);
    const Eigen::VectorXd invVarR = Eigen::VectorXd::LinSpaced(nobs, 1.0, nobs);
    const ObsErrorBogus R(nobs);
    Eigen::VectorXf wa(nana);
    Eigen::MatrixXf Wa(nana, nens);

    const auto times = oops::detGETKF_computeWeights(
        dy.cast<float>(), Yb, YbOrig, R, infl, svd, wa, Wa);

    // LAPACK implementation

    constexpr int getkf_inflation = 0;
    constexpr int denkf = 0;
    constexpr int getkf = 1;

    Eigen::VectorXf wa_f(nana);
    Eigen::MatrixXf Wa_f(nana, nens);

    const auto tL0 = std::chrono::system_clock::now();

    // cast eigen<double> to eigen<float>
    const Eigen::VectorXf dy_f = dy.cast<float>();
    const Eigen::MatrixXf Yb_f = Yb.cast<float>();
    const Eigen::MatrixXf YbOrig_f = YbOrig.cast<float>();
    const Eigen::VectorXf invVarR_f = invVarR.cast<float>();

    oops::letkf_core_f90(nobs,
                         Yb_f.data(),
                         YbOrig_f.data(),
                         dy_f.data(),
                         wa_f.data(),
                         Wa_f.data(),
                         invVarR_f.data(),
                         nana,
                         neig,
                         getkf_inflation,
                         denkf,
                         getkf,
                         infl);
    const auto tL1 = std::chrono::system_clock::now();

    // Timing information

    oops::Log::info() << std::endl;
    oops::Log::info() << "GETKF timing (ms)" << std::endl;
    oops::Log::info() << std::endl;

    oops::Log::info() << "Eigen:" << std::endl;
    oops::Log::info() << " compute Y^T R^-1: " << timeDifference(times[0], times[1]) << std::endl;
    oops::Log::info() << " compute Y^T R^-1 Y + I: " << timeDifference(times[1], times[2])
                      << std::endl;
    oops::Log::info() << " eigendecomposition: " << timeDifference(times[2], times[3]) << std::endl;
    oops::Log::info() << " compute Pa: " << timeDifference(times[3], times[4]) << std::endl;
    oops::Log::info() << " compute wa: " << timeDifference(times[4], times[5]) << std::endl;
    oops::Log::info() << " compute Wa: " << timeDifference(times[5], times[6]) << std::endl;
    oops::Log::info() << "Eigen total: " << timeDifference(times[0], times[6]) << std::endl;
    oops::Log::info() << std::endl;

    oops::Log::info() << "LAPACK total: " << timeDifference(tL0, tL1) << std::endl;
    oops::Log::info() << std::endl;

    // Compare the weights produced by the two implementations

    compareWeights(wa, wa_f, Wa, Wa_f);
  }

  void eigendecompositionWithSVD(const int nobs, const int nens, const int neig) {
    srand(1);

    // Use fake data for eigendecomposition (evd) and singular value decomposition (svd)

    // GETKF - deterministic
    const int nana = nens * neig;
    constexpr float infl = 1.0;
    const Eigen::VectorXd dy_Gevd = Eigen::VectorXd::Random(nobs);
    const Eigen::MatrixXf YbOrig_Gevd = Eigen::MatrixXf::Random(nens, nobs);
    const Eigen::MatrixXf Yb_Gevd = Eigen::MatrixXf::Random(nana, nobs);
    const Eigen::VectorXd invVarR_Gevd = Eigen::VectorXd::LinSpaced(nobs, 1.0, nobs);
    const ObsErrorBogus R(nobs);
    Eigen::VectorXf wa_Gevd(nana);
    Eigen::MatrixXf Wa_Gevd(nana, nens);
    const auto dy_Gsvd = dy_Gevd;
    const auto YbOrig_Gsvd = YbOrig_Gevd;
    const auto Yb_Gsvd = Yb_Gevd;
    const auto invVarR_Gsvd = invVarR_Gevd;
    Eigen::VectorXf wa_Gsvd(nana);
    Eigen::MatrixXf Wa_Gsvd(nana, nens);

    bool svd = false;
    const auto GETKF_times_EigenValue = oops::detGETKF_computeWeights(
        dy_Gevd.cast<float>(), Yb_Gevd, YbOrig_Gevd, R,
        infl, svd, wa_Gevd, Wa_Gevd);
    svd = true;
    const auto GETKF_times_SingularValue = oops::detGETKF_computeWeights(
        dy_Gsvd.cast<float>(), Yb_Gsvd, YbOrig_Gsvd, R,
        infl, svd, wa_Gsvd, Wa_Gsvd);

    // LETKF - deterministic
    const Eigen::VectorXd dy_Levd = Eigen::VectorXd::Random(nobs);
    const Eigen::MatrixXf Yb_Levd = Eigen::MatrixXf::Random(nens, nobs);
    const Eigen::VectorXd invVarR_Levd = Eigen::VectorXd::LinSpaced(nobs, 1.0, nobs);
    Eigen::VectorXf wa_Levd(nens);
    Eigen::MatrixXf Wa_Levd(nens, nens);
    const auto dy_Lsvd = dy_Levd;
    const auto Yb_Lsvd = Yb_Levd;
    const auto invVarR_Lsvd = invVarR_Levd;
    Eigen::VectorXf wa_Lsvd(nens);
    Eigen::MatrixXf Wa_Lsvd(nens, nens);
    svd = false;
    const auto LETKF_times_EigenValue = oops::detLETKF_computeWeights(
        dy_Levd.cast<float>(), Yb_Levd, R,
        (nens - 1) / 1.0, svd, wa_Levd, Wa_Levd);

    svd = true;
    const auto LETKF_times_SingularValue = oops::detLETKF_computeWeights(
        dy_Lsvd.cast<float>(), Yb_Lsvd, R,
        (nens - 1) / 1.0, svd, wa_Lsvd, Wa_Lsvd);

    // Timing information

    oops::Log::info() << std::endl;
    oops::Log::info() << "GETKF timing (ms)" << std::endl;
    oops::Log::info() << std::endl;

    oops::Log::info() << "Eigen value decomposition:" << std::endl;
    oops::Log::info() << " compute Y^T R^-1: "
                      << timeDifference(GETKF_times_EigenValue[0], GETKF_times_EigenValue[1])
                      << std::endl;
    oops::Log::info() << " compute Y^T R^-1 Y + I: "
                      << timeDifference(GETKF_times_EigenValue[1], GETKF_times_EigenValue[2])
                      << std::endl;
    oops::Log::info() << " eigendecomposition: "
                      << timeDifference(GETKF_times_EigenValue[2], GETKF_times_EigenValue[3])
                      << std::endl;
    oops::Log::info() << " compute Pa: "
                      << timeDifference(GETKF_times_EigenValue[3], GETKF_times_EigenValue[4])
                      << std::endl;
    oops::Log::info() << " compute wa: "
                      << timeDifference(GETKF_times_EigenValue[4], GETKF_times_EigenValue[5])
                      << std::endl;
    oops::Log::info() << " compute Wa: "
                      << timeDifference(GETKF_times_EigenValue[5], GETKF_times_EigenValue[6])
                      << std::endl;
    oops::Log::info() << "total: "
                      << timeDifference(GETKF_times_EigenValue[0], GETKF_times_EigenValue[6])
                      << std::endl;
    oops::Log::info() << std::endl;

    oops::Log::info() << "Singular value decomposition:" << std::endl;
    oops::Log::info() << " compute Y^T R^-1: "
                      << timeDifference(GETKF_times_SingularValue[0], GETKF_times_SingularValue[1])
                      << std::endl;
    oops::Log::info() << " compute Y^T R^-1 Y + I: "
                      << timeDifference(GETKF_times_SingularValue[1], GETKF_times_SingularValue[2])
                      << std::endl;
    oops::Log::info() << " SVD: "
                      << timeDifference(GETKF_times_SingularValue[2], GETKF_times_SingularValue[3])
                      << std::endl;
    oops::Log::info() << " compute Pa: "
                      << timeDifference(GETKF_times_SingularValue[3], GETKF_times_SingularValue[4])
                      << std::endl;
    oops::Log::info() << " compute wa: "
                      << timeDifference(GETKF_times_SingularValue[4], GETKF_times_SingularValue[5])
                      << std::endl;
    oops::Log::info() << " compute Wa: "
                      << timeDifference(GETKF_times_SingularValue[5], GETKF_times_SingularValue[6])
                      << std::endl;
    oops::Log::info() << "total: "
                      << timeDifference(GETKF_times_SingularValue[0], GETKF_times_SingularValue[6])
                      << std::endl;
    oops::Log::info() << std::endl;

    oops::Log::info() << std::endl;
    oops::Log::info() << "LETKF timing (ms)" << std::endl;
    oops::Log::info() << std::endl;

    oops::Log::info() << std::endl;
    oops::Log::info() << "Eigen value decomposition:" << std::endl;
    oops::Log::info() << " compute Y^T R^-1 Y: "
                      << timeDifference(LETKF_times_EigenValue[0], LETKF_times_EigenValue[1])
                      << std::endl;
    oops::Log::info() << " eigendecomposition: "
                      << timeDifference(LETKF_times_EigenValue[1], LETKF_times_EigenValue[2])
                      << std::endl;
    oops::Log::info() << " compute Pa: "
                      << timeDifference(LETKF_times_EigenValue[2], LETKF_times_EigenValue[3])
                      << std::endl;
    oops::Log::info() << " compute weights: "
                      << timeDifference(LETKF_times_EigenValue[3], LETKF_times_EigenValue[4])
                      << std::endl;
    oops::Log::info() << "total: "
                      << timeDifference(LETKF_times_EigenValue[0], LETKF_times_EigenValue[4])
                      << std::endl;
    oops::Log::info() << std::endl;

    oops::Log::info() << "Singular value decomposition:" << std::endl;
    oops::Log::info() << " compute Y^T R^-1 Y: "
                      << timeDifference(LETKF_times_SingularValue[0], LETKF_times_SingularValue[1])
                      << std::endl;
    oops::Log::info() << " SVD: "
                      << timeDifference(LETKF_times_SingularValue[1], LETKF_times_SingularValue[2])
                      << std::endl;
    oops::Log::info() << " compute Pa: "
                      << timeDifference(LETKF_times_SingularValue[2], LETKF_times_SingularValue[3])
                      << std::endl;
    oops::Log::info() << " compute weights: "
                      << timeDifference(LETKF_times_SingularValue[3], LETKF_times_SingularValue[4])
                      << std::endl;
    oops::Log::info() << "total: "
                      << timeDifference(LETKF_times_SingularValue[0], LETKF_times_SingularValue[4])
                      << std::endl;
    oops::Log::info() << std::endl;

    // Compare the weights produced by the two implementations
    compareWeights(wa_Lsvd, wa_Levd, Wa_Lsvd, Wa_Levd);
    compareWeights(wa_Gsvd, wa_Gevd, Wa_Gsvd, Wa_Gevd);
  }

  void test_ETKF(const eckit::LocalConfiguration & conf) {
    // Test array sizes

    const int nobs = conf.getInt("nobs");
    const int nens = conf.getInt("nens");
    const int neig = conf.getInt("neig");

    // LETKF test

    LETKF(nobs, nens);

    // GETKF test

    GETKF(nobs, nens, neig);

    // svd vs eigendecomposition test
    eigendecompositionWithSVD(nobs, nens, neig);
  }

  class ETKF : public oops::Test {
   public:
    using oops::Test::Test;
   private:
    std::string testid() const override {return "test::ETKF";}
    void register_tests() const override {
      std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

      const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
      for (const std::string & testCaseName : conf.keys())
        {
          const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(),
                                                       testCaseName);
          ts.emplace_back(CASE("ETKF/" + testCaseName, testCaseConf)
                          {
                            test_ETKF(testCaseConf);
                          });
        }
    }
    void clear() const override {}
  };

}  // namespace test

#endif  // TEST_ASSIMILATION_ETKF_H_
