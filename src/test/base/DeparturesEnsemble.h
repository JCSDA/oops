/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_BASE_DEPARTURESENSEMBLE_H_
#define TEST_BASE_DEPARTURESENSEMBLE_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename OBS> void testDeparturesEnsemble() {
  typedef ObsTestsFixture<OBS>            Test_;
  typedef oops::DeparturesEnsemble<OBS>   DeparturesEnsemble_;

  // tests depaturesEnsemle.packEigen()
  // test verifies that computing rms by-hand on the packed array agrees with
  // rms for each ensemble member computed using depatures.rms() method
  size_t myNens = 5;
  DeparturesEnsemble_ yens(Test_::obspace(), myNens);
  size_t myNobs = yens[0].nobs();

  std::vector<double> rms1(myNens);
  for (size_t ii=0; ii < myNens; ++ii) {
    yens[ii].random();
    rms1[ii] = yens[ii].rms();
  }

  Eigen::MatrixXd yEnsPack = yens.packEigen();
  Eigen::VectorXd rms2 = yEnsPack.rowwise().norm() / sqrt(myNobs);
  oops::Log::info() << "norm of eigen stuff: " << rms2 << std::endl;
  for (size_t iens = 0; iens < myNens; ++iens) {
    oops::Log::test() << "ii=" << iens << " yens[ii].rms=" << rms1[iens]
                      <<  " rms(ypack[ii,:])=" << rms2[iens] << std::endl;
    EXPECT(oops::is_close(rms2[iens], rms1[iens], 1.e-14));
  }
}

template <typename OBS> class DeparturesEnsemble : public oops::Test {
 public:
  DeparturesEnsemble() {}
  virtual ~DeparturesEnsemble() {}
 private:
  std::string testid() const {return "test::DeparturesEnsemble<" + OBS::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();
    ts.emplace_back(CASE("base/DeparturesEnsemble/testDeparturesEnsemble")
      { testDeparturesEnsemble<OBS>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_BASE_DEPARTURESENSEMBLE_H_
