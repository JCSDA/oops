/*
 * (C) Crown copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_ASSIMILATION_SUBENSEMBLESPLITTER_H_
#define TEST_ASSIMILATION_SUBENSEMBLESPLITTER_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <string>
#include <tuple>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/assimilation/SubensembleSplitter.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"

namespace test {

// -----------------------------------------------------------------------------
void testKGO() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  std::vector<eckit::LocalConfiguration> testCases;
  conf.get("test cases", testCases);
  const std::vector<std::string> validCases = {"Contiguous",
                                               "Random, no seed",
                                               "Random, seeded",
                                               "Contiguous (modulated)",
                                               "Random, no seed (modulated)",
                                               "Random, seeded (modulated)"};

  const size_t nens = 6;
  for (const eckit::LocalConfiguration & testCase : testCases) {
    eckit::LocalConfiguration testConf;
    testCase.get("test case", testConf);
    const std::string testName = testConf.getString("name");
    if (!(std::find(validCases.begin(), validCases.end(), testName) != validCases.end())) {
      continue;
    }
    oops::Log::info() << "Running test case: " << testName << "\n";

    // Setup splitter and split
    eckit::LocalConfiguration splitterConf = testConf.getSubConfiguration("subensemble splitter");
    const size_t neig = testConf.getUnsigned("number of eigenvalues", 1);
    oops::SubensembleSplitter splitter(nens, neig, splitterConf);
    splitter.split();

    // Get the exclusion components
    const size_t excludedSubensemble = testConf.getUnsigned("excluded subensemble", 0);
    const bool modulated = (neig != 1);
    std::tuple<Eigen::SparseMatrix<float>, Eigen::SparseMatrix<float>, std::vector<size_t>>
    projectionMatrices = splitter.getProjectionMatrices(excludedSubensemble, modulated);
    const Eigen::SparseMatrix<float> & excludedProjection = std::get<0>(projectionMatrices);
    const std::vector<size_t> & excludedMembers = std::get<2>(projectionMatrices);
    for (size_t member : excludedMembers) {
      oops::Log::info() << "Member = " << member << " excluded" << "\n";
    }

    // Testing KGO
    const size_t nsubens = splitterConf.getUnsigned("number of subensembles");
    const size_t nhat = nens*neig - (nens/nsubens)*neig;
    const std::vector<size_t> KGOLocs = testConf.getUnsignedVector("location of ones");
    const std::vector<size_t> KGOMems = testConf.getUnsignedVector("excluded members");
    float sumOfElem = 0.0f;
    const float tol = 1e-10f;
    const size_t nrows = excludedProjection.rows();
    const size_t ncols = excludedProjection.cols();
    for (size_t col = 0; col < ncols; ++col) {
      const size_t locRow = KGOLocs[col];
      EXPECT(oops::is_close(excludedProjection.coeff(locRow, col), 1.0f, tol));

      for (size_t row = 0; row < nrows; ++row) {
        sumOfElem += excludedProjection.coeff(row, col);
      }
    }
    EXPECT(oops::is_close(sumOfElem, static_cast<float>(nhat), tol));
    EXPECT(std::is_permutation(excludedMembers.begin(), excludedMembers.end(),
                               KGOMems.begin(), KGOMems.end()));
  }
}

// -----------------------------------------------------------------------------
void testErrorHandling() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  std::vector<eckit::LocalConfiguration> testCases;
  conf.get("test cases", testCases);
  const std::vector<std::string> validCases = {"Invalid number of eigenvalues",
                                               "Too few subensembles",
                                               "Too many subensembles",
                                               "Does not divide",
                                               "Excluded subensemble too large",
                                               "Invalid splitting method"};

  const size_t nens = 6;
  for (const eckit::LocalConfiguration & testCase : testCases) {
    eckit::LocalConfiguration testConf;
    testCase.get("test case", testConf);
    const std::string testName = testConf.getString("name");
    if (!(std::find(validCases.begin(), validCases.end(), testName) != validCases.end())) {
      continue;
    }
    oops::Log::info() << "Running test case: " << testName << "\n";

    // Testing error handling
    eckit::LocalConfiguration splitterConf = testConf.getSubConfiguration("subensemble splitter");
    const std::string msg = testConf.getString("expectExceptionWithMessage");
    if ((testName == "Invalid splitting method")) {
      oops::SubensembleSplitter splitter(nens, splitterConf);
      EXPECT_THROWS_MSG(splitter.split(), msg.c_str());
    } else if ((testName == "Excluded subensemble too large")) {
      oops::SubensembleSplitter splitter(nens, splitterConf);
      splitter.split();
      size_t excludedSubensemble = testConf.getUnsigned("excluded subensemble");
      EXPECT_THROWS_MSG(splitter.getProjectionMatrices(excludedSubensemble, false), msg.c_str());
    } else {
      const size_t neig = testConf.getUnsigned("number of eigenvalues", 1);
      EXPECT_THROWS_MSG(oops::SubensembleSplitter splitter(nens, neig, splitterConf),
                        msg.c_str());
    }
  }
}

// -----------------------------------------------------------------------------
void testApplyProjection() {
  const eckit::LocalConfiguration conf(::test::TestEnvironment::config());
  std::vector<eckit::LocalConfiguration> testCases;
  conf.get("test cases", testCases);
  const std::vector<std::string> validCases = {"Unmodulated projection demo",
                                               "Modulated projection demo",
                                               "Unmodulated projection demo random",
                                               "Modulated projection demo random"};

  const size_t nens = 3;
  for (const eckit::LocalConfiguration & testCase : testCases) {
    eckit::LocalConfiguration testConf;
    testCase.get("test case", testConf);
    const std::string testName = testConf.getString("name");
    if (!(std::find(validCases.begin(), validCases.end(), testName) != validCases.end())) {
      continue;
    }
    oops::Log::info() << "Running test case: " << testName << "\n";

    // Demonstrating application of projection matrices
    const size_t neig = testConf.getUnsigned("number of eigenvalues", 1);
    eckit::LocalConfiguration splitterConf = testConf.getSubConfiguration("subensemble splitter");
    const size_t nsubens = splitterConf.getUnsigned("number of subensembles");
    const size_t nanal = nens*neig;
    const size_t nhat = nanal - (nens/nsubens)*neig;
    const Eigen::MatrixXf mat1 = Eigen::MatrixXf::Random(nanal, nanal);
    const Eigen::MatrixXf mat2 = Eigen::MatrixXf::Random(nens, nens);

    // Setup splitter and split
    oops::SubensembleSplitter splitter(nens, neig, splitterConf);
    splitter.split();

    // Get the exclusion components
    std::tuple<Eigen::SparseMatrix<float>, Eigen::SparseMatrix<float>, std::vector<size_t>>
    projectionMatrices;
    const bool modulated = (neig != 1);
    Eigen::MatrixXf sumOfComplements = Eigen::MatrixXf::Zero(nens, nens);
    const float tol = 1e-10;
    for (size_t isubens = 0; isubens < nsubens; ++isubens) {
      projectionMatrices = splitter.getProjectionMatrices(isubens, modulated);
      const Eigen::SparseMatrix<float> & excludedProjection = std::get<0>(projectionMatrices);
      const Eigen::SparseMatrix<float> & includedProjection = std::get<1>(projectionMatrices);
      const std::vector<size_t> & excludedMembers = std::get<2>(projectionMatrices);

      // Apply projection
      // Project out the current subensemble members
      const Eigen::MatrixXf res1 = mat1*excludedProjection;
      // Apply update to current ensemble members
      const Eigen::MatrixXf res2 = mat2*includedProjection;

      // Test results
      /* Our matrix result for res1 should consist of the columns
         not excluded in original matrix */
      const size_t nrows1 = res1.rows();
      const size_t ncols1 = res1.cols();
      EXPECT(nrows1 == nanal);
      EXPECT(ncols1 == nhat);
      size_t resCol = 0;
      for (size_t iens = 0; iens < nens; ++iens) {
        const bool isExcluded = (std::find(excludedMembers.begin(),
                                           excludedMembers.end(),
                                           iens) != excludedMembers.end());
        if (!isExcluded) {
          for (size_t ieig = 0; ieig < neig; ++ieig) {
            const size_t origCol = iens*neig + ieig;
            for (size_t row = 0; row < nanal; ++row) {
              EXPECT(oops::is_close(mat1(row, origCol),
                                    res1(row, resCol),
                                    tol));
            }
            ++resCol;
          }
        }
      }

      /* Our matrix result for res2 should consist of the columns
         excluded in original matrix, and zeroes elsewhere */
      const size_t nrows2 = res2.rows();
      const size_t ncols2 = res2.cols();
      EXPECT(nrows2 == nens);
      EXPECT(ncols2 == nens);
      for (size_t col = 0; col < nens; ++col) {
        const bool isExcluded = (std::find(excludedMembers.begin(),
                                           excludedMembers.end(),
                                           col) != excludedMembers.end());
        for (size_t row = 0; row < nens; ++row) {
          const float kgo = (isExcluded) ? mat2(row, col) : 0.0;
          EXPECT(oops::is_close(kgo,
                                res2(row, col),
                                tol));
        }
      }
      sumOfComplements += includedProjection;
    }

    /* Our matrix result for the sum of complements should be identity,
       in other words no subensemble contains a member present in the others */
    const Eigen::MatrixXf sumKGO = Eigen::MatrixXf::Identity(nens, nens);
    for (size_t row = 0; row < nens; ++row) {
      for (size_t col = 0; col < nens; ++col) {
        EXPECT(oops::is_close(sumKGO(row, col),
                              sumOfComplements(row, col),
                              tol));
      }
    }
  }
}

// -----------------------------------------------------------------------------
class SubensembleSplitter : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::SubensembleSplitter";}
  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("oops/SubensembleSplitter/")
                    {
                      testKGO();
                    });

    ts.emplace_back(CASE("oops/SubensembleSplitter/")
                    {
                      testErrorHandling();
                    });

    ts.emplace_back(CASE("oops/SubensembleSplitter/")
                    {
                      testApplyProjection();
                    });
  }
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_ASSIMILATION_SUBENSEMBLESPLITTER_H_
