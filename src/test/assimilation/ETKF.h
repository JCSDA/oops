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

  void compareWeights(const Eigen::VectorXd & wa, const Eigen::VectorXd & wa_d,
                      const Eigen::MatrixXd & Wa, const Eigen::MatrixXd & Wa_d) {
    for (int i = 0; i < wa.rows(); ++i) {
      EXPECT(oops::is_close_absolute(wa(i), wa_d(i), 1.0e-6));
    }
    for (int i = 0; i < Wa.rows(); ++i) {
      for (int j = 0; j < Wa.cols(); ++j) {
        EXPECT(oops::is_close_absolute(Wa(i, j), Wa_d(i, j), 1.0e-5));
      }
    }
  }

  void LETKF(const int nobs, const int nens) {
    srand(1);

    const Eigen::VectorXd dy = Eigen::VectorXd::Random(nobs);
    const Eigen::MatrixXf Yb = Eigen::MatrixXf::Random(nens, nobs);
    const Eigen::VectorXd invVarR = Eigen::VectorXd::LinSpaced(nobs, 1.0, nobs);
    Eigen::VectorXd wa(nens);
    Eigen::MatrixXd Wa(nens, nens);

    // Eigen implementation

    const auto times = oops::detLETKF_computeWeights(dy, Yb, invVarR, (nens - 1) / 1.0, wa, Wa);

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

    const Eigen::VectorXd wa_d = wa_f.cast<double>();
    const Eigen::MatrixXd Wa_d = Wa_f.cast<double>();

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

    compareWeights(wa, wa_d, Wa, Wa_d);
  }

  void GETKF(const int nobs, const int nens, const int neig) {
    srand(1);

    const int nana = nens * neig;
    constexpr float infl = 1.0;
    const Eigen::VectorXd dy = Eigen::VectorXd::Random(nobs);
    const Eigen::MatrixXf YbOrig = Eigen::MatrixXf::Random(nens, nobs);
    // modulated ensemble
    const Eigen::MatrixXf Yb = Eigen::MatrixXf::Random(nana, nobs);
    const Eigen::VectorXd invVarR = Eigen::VectorXd::LinSpaced(nobs, 1.0, nobs);
    Eigen::VectorXd wa(nana);
    Eigen::MatrixXd Wa(nana, nens);

    const auto times = oops::detGETKF_computeWeights(dy, Yb, YbOrig, invVarR, infl, wa, Wa);

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

    const Eigen::VectorXd wa_d = wa_f.cast<double>();
    const Eigen::MatrixXd Wa_d = Wa_f.cast<double>();

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

    compareWeights(wa, wa_d, Wa, Wa_d);
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
