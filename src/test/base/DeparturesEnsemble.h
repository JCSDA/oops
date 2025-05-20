/*
 * (C) Copyright 2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_BASE_DEPARTURESENSEMBLE_H_
#define TEST_BASE_DEPARTURESENSEMBLE_H_

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/testing/Test.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/ObsSpaces.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {
/// \brief Tests DeparturesEnsemble class.
/// \details
/// 1. Create DeparturesEnsemble object (depsEns) by sequentially defining
///    nens (ensemble size) Departures objects (deps)
/// 2. Calculate and compare rms of depsEns and nens deps objects.
/// 3. depsRms and depsEnsRms needs to be identical to pass this test.

template <typename OBS> void testDeparturesEnsemble() {
  typedef ObsTestsFixture<OBS>                   Test_;
  typedef oops::Departures<OBS>                  Departures_;
  typedef oops::DeparturesEnsemble<OBS>          DeparturesEnsemble_;
  typedef oops::ObsSpaces<OBS>                   ObsSpaces_;

  // configure obs spaces
  const eckit::LocalConfiguration obsConfig(TestEnvironment::config(), "observations");
  const ObsSpaces_ & obspaces = Test_::obspace();

  // declare nens, Deps, DepsEns
  const size_t nens = 4;
  Departures_ deps(obspaces);
  DeparturesEnsemble_ depsEns(obspaces, nens);

  // define mask departure object and set its default values to zero
  // it will be used in the packEigen method.
  Departures_ mask(obspaces);
  mask.zero();
  const Eigen::VectorXd dep0 = deps.packEigen(mask);
  const size_t nobs = dep0.size();

  if ( nobs < 1 ) {
    oops::Log::error() << "nobs cannot be less than 1";
    throw eckit::BadParameter("nobs cannot be less than 1");
  }

  // initialize nens deps randomly to generate depsEns
  // and calcuate rms and apply packEigen to deps
  Eigen::MatrixXd depsMat(nens, nobs);
  Eigen::MatrixXf depsEnsMat(nens, nobs);
  std::vector<float> depsRms(nens);
  std::vector<float> depsEnsRms(nens);

  for (size_t iens =  0; iens < nens; ++iens) {
    deps.random();
    depsMat.row(iens) = deps.packEigen(mask);
    depsRms[iens] = static_cast<float>(deps.rms());
    depsEns.setData(iens, deps);
  }

  // apply packEigen to depsEns and calculate rms
  depsEnsMat = depsEns.packEigen(mask);
  for (size_t iens =  0; iens < nens; ++iens) {
    depsEnsRms[iens] = static_cast<float>(depsEns.getData(iens).rms());
  }

  // compare rms vector and packEigen matrix in nens deps
  // and depsEns
  for (size_t i = 0; i < nens; i++) {
    oops::Log::debug() << "Departures and DeparturesEnsembles RMS at mem " << i << ":" << std::endl;
    oops::Log::debug() << depsRms[i] << ", " << depsEnsRms[i] << std::endl;
    EXPECT_EQUAL(depsRms[i], depsEnsRms[i]);
  }

  // compare depsMat and depsEnsMat to check DeparturesEnsemble::packEigen()
  const float tol = 1e-6f;
  for (size_t iens = 0; iens < nens; iens ++) {
    EXPECT(oops::is_close(static_cast<float>(depsMat.row(iens).squaredNorm()),
                          depsEnsMat.row(iens).squaredNorm(), tol));
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> class DeparturesEnsemble : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;

 public:
  using oops::Test::Test;
  virtual ~DeparturesEnsemble() = default;
 private:
  std::string testid() const override {return "test::DeparturesEnsemble<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("base/DeparturesEnsemble/testDeparturesEnsemble")
      { testDeparturesEnsemble<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_BASE_DEPARTURESENSEMBLE_H_
