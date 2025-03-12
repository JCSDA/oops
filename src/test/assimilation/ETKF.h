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
#include "oops/assimilation/gletkfInterface.h"
#include "oops/assimilation/LETKFSolver.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"

namespace test {

  auto timeDifference = [](auto tA, auto tB) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(tB - tA).count();
  };

  void compareWeights(const Eigen::VectorXf & wa, const Eigen::VectorXf & wa_f,
                      const Eigen::MatrixXf & Wa, const Eigen::MatrixXf & Wa_f) {
    for (int i = 0; i < wa.rows(); ++i) {
      EXPECT(oops::is_close_absolute(wa(i), wa_f(i), 1.0e-6f));
    }
    for (int i = 0; i < Wa.rows(); ++i) {
      for (int j = 0; j < Wa.cols(); ++j) {
        EXPECT(oops::is_close_absolute(Wa(i, j), Wa_f(i, j), 1.0e-5f));
      }
    }
  }

  void LETKF(const int nobs, const int nens) {
    srand(1);

    const Eigen::VectorXf dy_f = Eigen::VectorXf::Random(nobs);
    const Eigen::MatrixXf Yb_f = Eigen::MatrixXf::Random(nens, nobs);
    const Eigen::VectorXf diagInvR_f = Eigen::VectorXf::LinSpaced(nobs, 1.0, nobs);

    // Eigen implementation

    const auto tE0 = std::chrono::system_clock::now();

    // work = Y^T R^-1 Y + (nens-1)/infl I
    Eigen::MatrixXf work = Yb_f * (diagInvR_f.asDiagonal() * Yb_f.transpose());
    work.diagonal() += Eigen::VectorXf::Constant(nens, nens - 1);

    const auto tE1 = std::chrono::system_clock::now();

    // eigenvalues and eigenvectors of the above matrix
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(work);
    const Eigen::VectorXf eival = es.eigenvalues().real();
    const Eigen::MatrixXf eivec = es.eigenvectors().real();

    const auto tE2 = std::chrono::system_clock::now();

    // Pa = [ Yb^T R^-1 Yb + (nens-1)/infl I ] ^-1
    work = eivec * eival.cwiseInverse().asDiagonal() * eivec.transpose();

    const auto tE3 = std::chrono::system_clock::now();

    // Wa = sqrt[(nens-1) Pa]
    const Eigen::MatrixXf Wa =
      eivec *
      ((nens - 1) * eival.array().inverse()).sqrt().matrix().asDiagonal() *
      eivec.transpose();

    // wa = Pa Yb^T R^-1 dy
    const Eigen::VectorXf wa = work * (Yb_f * (diagInvR_f.asDiagonal() * dy_f));

    const auto tE4 = std::chrono::system_clock::now();

    // LAPACK implementation

    constexpr int neig = 1;
    constexpr int getkf_inflation = 0;
    constexpr int denkf = 0;
    constexpr int getkf = 0;
    constexpr float infl = 1.0;

    Eigen::VectorXf wa_f(nens);
    Eigen::MatrixXf Wa_f(nens, nens);

    const auto tL0 = std::chrono::system_clock::now();

    oops::letkf_core_f90(nobs,
                         Yb_f.data(),
                         Yb_f.data(),
                         dy_f.data(),
                         wa_f.data(),
                         Wa_f.data(),
                         diagInvR_f.data(),
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
    oops::Log::info() << " compute Y^T R^-1 Y: " << timeDifference(tE0, tE1) << std::endl;
    oops::Log::info() << " eigendecomposition: " << timeDifference(tE1, tE2) << std::endl;
    oops::Log::info() << " compute Pa: " << timeDifference(tE2, tE3) << std::endl;
    oops::Log::info() << " compute weights: " << timeDifference(tE3, tE4) << std::endl;
    oops::Log::info() << "Eigen total: " << timeDifference(tE0, tE4) << std::endl;
    oops::Log::info() << std::endl;
    oops::Log::info() << "LAPACK total: " << timeDifference(tL0, tL1) << std::endl;
    oops::Log::info() << std::endl;

    // Compare the weights produced by the two implementations

    compareWeights(wa, wa_f, Wa, Wa_f);
  }

  void GETKF(const int nobs, const int nens, const int neig) {
    srand(1);

    const int nana = nens * neig;
    constexpr float infl = 1.0;
    const Eigen::VectorXf dy_f = Eigen::VectorXf::Random(nobs);
    const Eigen::MatrixXf YbOrig_f = Eigen::MatrixXf::Random(nens, nobs);
    // modulated ensemble
    const Eigen::MatrixXf Yb_f = Eigen::MatrixXf::Random(nana, nobs);
    const Eigen::VectorXf diagInvR_f = Eigen::VectorXf::LinSpaced(nobs, 1.0, nobs);

    // Eigen implementation

    const auto tE0 = std::chrono::system_clock::now();

    // Identity
    const Eigen::VectorXf I = Eigen::VectorXf::Constant(nana, 1.0);

    // Yb R^-1
    const Eigen::MatrixXf YbRinv = Yb_f * diagInvR_f.asDiagonal();

    const auto tE1 = std::chrono::system_clock::now();

    // Yb R^-1 Yb^T + (nens - 1) I / infl
    const Eigen::MatrixXf YbRinvYbpI =
      YbRinv * Yb_f.transpose() +
      I.asDiagonal().toDenseMatrix() * (nens - 1) / infl;

    const auto tE2 = std::chrono::system_clock::now();

    // Eigendecomposition
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(YbRinvYbpI);
    const Eigen::VectorXf eival = es.eigenvalues().real();
    const Eigen::MatrixXf eivec = es.eigenvectors().real();

    const auto tE3 = std::chrono::system_clock::now();

    // Pa = [Yb R^-1 Yb^T + (nens - 1)/infl I]^-1
    const Eigen::MatrixXf Pa =
      eivec * eival.cwiseInverse().asDiagonal() * eivec.transpose();

    const auto tE4 = std::chrono::system_clock::now();

    // wa_f = Pa Yb R^-1 dy
    const Eigen::VectorXf wa = Pa * (YbRinv * dy_f);

    const auto tE5 = std::chrono::system_clock::now();

    // Normalisation
    const float norm = 1.0 / (nens - 1.0);

    // (I - Gamma^{-1/2} / (nens - 1)) * (Gamma - (nens - 1) I / rho)
    const Eigen::VectorXf diag = (I - (norm * eival).cwiseInverse().cwiseSqrt()).
      cwiseProduct((eival - I / (infl * norm)).cwiseInverse());

    // C ((I - Gamma^{-1/2} / (nens - 1)) * (Gamma - (nens - 1) I / rho)) C^T
    const Eigen::MatrixXf scaledCCT = eivec * diag.asDiagonal() * eivec.transpose();

    const auto tE6 = std::chrono::system_clock::now();

    // Yb R^-1 YbOrig^T
    const Eigen::MatrixXf YbRinvYbOrig = YbRinv * YbOrig_f.transpose();

    const auto tE7 = std::chrono::system_clock::now();

    // Wa_f
    const Eigen::MatrixXf Wa = -scaledCCT * YbRinvYbOrig;

    const auto tE8 = std::chrono::system_clock::now();

    // LAPACK implementation

    constexpr int getkf_inflation = 0;
    constexpr int denkf = 0;
    constexpr int getkf = 1;

    Eigen::VectorXf wa_f(nana);
    Eigen::MatrixXf Wa_f(nana, nens);

    const auto tL0 = std::chrono::system_clock::now();

    oops::letkf_core_f90(nobs,
                         Yb_f.data(),
                         YbOrig_f.data(),
                         dy_f.data(),
                         wa_f.data(),
                         Wa_f.data(),
                         diagInvR_f.data(),
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
    oops::Log::info() << " compute Y^T R^-1: " << timeDifference(tE0, tE1) << std::endl;
    oops::Log::info() << " compute Y^T R^-1 Y + I: " << timeDifference(tE1, tE2) << std::endl;
    oops::Log::info() << " eigendecomposition: " << timeDifference(tE2, tE3) << std::endl;
    oops::Log::info() << " compute Pa: " << timeDifference(tE3, tE4) << std::endl;
    oops::Log::info() << " compute wa: " << timeDifference(tE4, tE5) << std::endl;
    oops::Log::info() << " compute scaled C C^T: " << timeDifference(tE5, tE6) << std::endl;
    oops::Log::info() << " compute Yb R^-1 YbOrig: " << timeDifference(tE6, tE7) << std::endl;
    oops::Log::info() << " compute Wa: " << timeDifference(tE7, tE8) << std::endl;
    oops::Log::info() << "Eigen total: " << timeDifference(tE0, tE8) << std::endl;
    oops::Log::info() << std::endl;

    oops::Log::info() << "LAPACK total: " << timeDifference(tL0, tL1) << std::endl;
    oops::Log::info() << std::endl;

    // Compare the weights produced by the two implementations

    compareWeights(wa, wa_f, Wa, Wa_f);
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
