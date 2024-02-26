/*
 * (C) Crown copyright 2021, Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

/// \file ObsErrorDiagZeroMeanPerturbations.cc
///
/// Verifies the 'zero-mean perturbations' option of ObsErrorDiag works correctly, i.e.
/// perturbations generated with this option have zero ensemble mean and the expected variance.

#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/ObsTable.h"
#include "oops/runs/Run.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/sqr.h"
#include "oops/util/TimeWindow.h"
#include "test/TestEnvironment.h"

namespace test {

CASE("test_obserrordiag_zeromeanpert") {
  std::vector<std::string> files = TestEnvironment::config().getStringVector("files");

  std::vector<double> originalObs;
  std::vector<double> obsErrors;
  std::vector<double> perturbedObs;
  std::vector<double> perturbedObsMean;

  // sum of ((perturbedObs[i] - originalObs[i]) / obsErrors[i])^2 over all locations i
  // in all ensemble members
  double normalizedVariance = 0.0;
  // nominal amplitude of perturbations; we assume it's identical for all members
  double randomAmplitude = 1.0;  // 1.0 -- default value

  size_t member = 0;
  for (const std::string &file : files) {
    ++member;

    const eckit::PathName confPath = file;
    const eckit::YAMLConfiguration memberConf(confPath);
    const eckit::LocalConfiguration costFunctionConf(memberConf, "cost function");
    const std::vector<eckit::LocalConfiguration> obsConfs =
        costFunctionConf.getSubConfigurations("observations.observers");
    const eckit::LocalConfiguration timeWindowConf(costFunctionConf, "time window");
    EXPECT_EQUAL(obsConfs.size(), 1);

    const eckit::LocalConfiguration obsErrorConf(obsConfs[0], "obs error");
    randomAmplitude = obsErrorConf.getDouble("obs perturbations amplitude");

    eckit::LocalConfiguration obsSpaceConf;
    obsSpaceConf.set("obsdatain.obsfile", obsConfs[0].getString("obs space.obsdataout.obsfile"));

    lorenz95::ObsTable obsSpace(obsSpaceConf, oops::mpi::world(),
                                util::TimeWindow(timeWindowConf),
                                oops::mpi::myself());
    obsSpace.getdb("ObsValue", originalObs);
    obsSpace.getdb("ObsError", obsErrors);
    obsSpace.getdb("EffectiveObsValue", perturbedObs);

    perturbedObsMean.resize(perturbedObs.size());
    for (size_t loc = 0; loc < perturbedObs.size(); ++loc)
      perturbedObsMean[loc] += perturbedObs[loc];

    for (size_t loc = 0; loc < perturbedObs.size(); ++loc)
      normalizedVariance += util::sqr((perturbedObs[loc] - originalObs[loc]) / obsErrors[loc]);
  }

  const size_t numMembers = member;

  // Check that the ensemble mean of perturbations is zero
  const double invNumMembers = 1.0 / files.size();
  for (size_t loc = 0; loc < perturbedObs.size(); ++loc)
    perturbedObsMean[loc] *= invNumMembers;

  // The input file contains numbers printed with 5 digits after the decimal point
  EXPECT(oops::are_all_close_absolute(perturbedObsMean, originalObs, 1e-5));

  // Check that perturbations have the expected variance
  normalizedVariance /= (perturbedObs.size() * numMembers);
  const double normalizedStdDev = std::sqrt(normalizedVariance);

  oops::Log::test() << "Std dev of perturbations / obs error: " << normalizedStdDev << "\n";
  oops::Log::test() << "Nominal value: " << randomAmplitude << std::endl;

  // If this test passes with this threshold, we can be reasonably sure the scaling factor
  // sqrt(n/(n-1)) included in zero-mean perturbations is correct
  double maxAcceptableRelativeError = 0.25 * std::sqrt(numMembers / (numMembers - 1));
  EXPECT(oops::is_close_relative(normalizedStdDev, randomAmplitude, maxAcceptableRelativeError));
}

class ObsErrorDiagZeroMeanPerturbations : public oops::Test {
 public:
  ObsErrorDiagZeroMeanPerturbations() {}

 private:
  std::string testid() const override {return "test::ObsErrorDiagZeroMeanPerturbations";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test

int main(int argc, char **argv) {
  oops::Run run(argc, argv);
  test::ObsErrorDiagZeroMeanPerturbations tests;
  return run.execute(tests);
}
