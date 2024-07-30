/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_ASSIMILATION_SPECTRALLMP_H_
#define TEST_ASSIMILATION_SPECTRALLMP_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/assimilation/SpectralLMP.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"

#include "test/assimilation/Vector3D.h"

namespace test {

  void test_SpectralLMP(const eckit::LocalConfiguration &conf)
  {
    oops::SpectralLMP<Vector3D, Vector3D> spectralLMP(conf);

    // assign vectors following DRPLanczosMinimizer.h
    std::vector<std::unique_ptr<Vector3D>> hvecs;
    std::vector<std::unique_ptr<Vector3D>> vvecs;
    std::vector<std::unique_ptr<Vector3D>> zvecs;
    std::vector<double> alphas;
    std::vector<double> betas;

    // Simple case
    hvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 1, 1)));
    hvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 1, 1)));
    hvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 1, 1)));
    vvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 1, 1)));
    vvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 1, 1)));
    vvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 1, 1)));
    zvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 1, 1)));
    zvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 1, 1)));
    zvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 1, 1)));

    alphas.push_back(1.0);
    alphas.push_back(1.0);
    betas.push_back(1.0);
    betas.push_back(1.0);

    spectralLMP.update(vvecs, hvecs, zvecs, alphas, betas);

    Vector3D pr(1, 1, 1);
    Vector3D zz(1, 1, 1);

    spectralLMP.updateObsBias(std::unique_ptr<Vector3D>(new Vector3D(1, 1, 1)));
    spectralLMP.multiply(pr, zz);

    // More complicated case
    hvecs.clear();
    vvecs.clear();
    zvecs.clear();
    alphas.clear();
    betas.clear();
    hvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(10, 1, -10)));
    hvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 10, 100)));
    hvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(-5, 5, 20)));
    vvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 3, 11)));
    vvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(-1, 2, 100)));
    vvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(6, 77, 7)));
    zvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 10, 1)));
    zvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(1, 21, 21)));
    zvecs.emplace_back(std::unique_ptr<Vector3D>(new Vector3D(100, 3, 70)));

    alphas.push_back(1000.0);
    alphas.push_back(33.0);
    betas.push_back(93.0);
    betas.push_back(2.0);

    spectralLMP.update(vvecs, hvecs, zvecs, alphas, betas);

    Vector3D pr2(1, 3, 1);
    Vector3D zz2(50, 1, 1);

    spectralLMP.multiply(pr2, zz2);
  }

  void test_SpectralLMP_Cmat(const eckit::LocalConfiguration &conf)
  {
    oops::SpectralLMP<Vector3D, Vector3D> spectralLMP(conf);

    Vector3D pr(1, 1, 1);
    Vector3D zz(1, 1, 1);

    EXPECT_THROWS_MSG(spectralLMP.multiply(pr, zz),
                      "The VarBC preconditioner matrix is not defined");
  }

  class SpectralLMP : public oops::Test {
   public:
    using oops::Test::Test;

   private:
    std::string testid() const override {return "test::SpectralLMP";}

    void register_tests() const override {
      std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

      const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
      for (const std::string & testCaseName : conf.keys())
        {
          const eckit::LocalConfiguration testCaseConf(::test::TestEnvironment::config(),
                                                       testCaseName);
          ts.emplace_back(CASE("SpectralLMP/" + testCaseName, testCaseConf)
                          {
                            test_SpectralLMP(testCaseConf);
                          });
          ts.emplace_back(CASE("SpectralLMPCmat/" + testCaseName, testCaseConf)
                          {
                            test_SpectralLMP_Cmat(testCaseConf);
                          });
        }
    }

    void clear() const override {}
  };

}  // namespace test

#endif  // TEST_ASSIMILATION_SPECTRALLMP_H_
