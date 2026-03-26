/*
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_GENERIC_VERTICALLOCEV_H_
#define TEST_GENERIC_VERTICALLOCEV_H_

#include <Eigen/Dense>

#include <cfloat>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/Variables.h"
#include "oops/generic/VerticalLocEV.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "test/interface/Increment.h"
#include "test/TestEnvironment.h"

namespace test {

// =============================================================================

template <typename MODEL> void testVerticalLocEV() {
  typedef IncrementFixture<MODEL>         Test_;
  typedef oops::Geometry<MODEL>           Geometry_;
  typedef oops::VerticalLocEV<MODEL>      VerticalLocEV_;
  typedef oops::Increment4D<MODEL>        Increment4D_;
  typedef oops::IncrementSet<MODEL>       IncrementSet_;
  typedef oops::State<MODEL>              State_;

  const Geometry_ & geometry = Test_::resol();
  eckit::LocalConfiguration vertlocconf(TestEnvironment::config(), "vertical localization");

  // make an empty state vector to be used to intitialize VerticalLocEV_
  State_ x(Test_::resol(), Test_::ctlvars(), Test_::time());

  VerticalLocEV_ vertloc(vertlocconf, x, x.variables());
  oops::Log::info() << "Number of eigenvalues used in VerticalLoc: " << vertloc.neig() << std::endl;

  //--- check for expected number of eigen modes
  int nEigExpected = TestEnvironment::config().getInt("expected neig");
  int neig = vertloc.neig();
  EXPECT(nEigExpected == neig);

  //--- check that truncation and rescaling was done correctly
  EXPECT(vertloc.testTruncateEvecs(geometry));

  //--- check modulation
  // if the increment dx=1, then modulated ensemble will be the same as (scaled) eigen vectors
  // then the dot products will satisfy orthogonality condition

  // need at least 2 eigen vectors for the following check to work.
  EXPECT(neig > 1);

  // construct instances of Increment4D_ and IncrementEnsemble_
  std::vector<util::DateTime> times;
  times.push_back(Test_::time());
  Increment4D_ dx1(Test_::resol(), Test_::ctlvars(), times);
  Increment4D_ dx2(Test_::resol(), Test_::ctlvars(), times);
  // Build an IncrementSet with one time and neig members
  std::vector<int> members(neig);
  std::iota(members.begin(), members.end(), 0);
  IncrementSet_ incEns(Test_::resol(), Test_::ctlvars(), times, oops::mpi::myself(), members);
  // set incEns to zero
  for (int i = 1; i < neig; ++i) {incEns(0, i).zero();}
  for (int i = 1; i < neig; ++i) {
    double n = incEns(0, 0).dot_product_with(incEns(0, i));
    EXPECT(n < 100*DBL_EPSILON);
  }

  // make increment of ones and check that it is indeed full of ones
  dx1.ones();
  dx2.random();
  double normRandom = dx2.dot_product_with(dx2);
  dx2.schur_product_with(dx1);
  EXPECT(std::abs(dx2.dot_product_with(dx2)-normRandom) < normRandom*DBL_EPSILON);

  // modulate increments
  vertloc.modulateIncrement(dx1, incEns);

  // check the orthogonality condition
  double n0 = incEns(0, 0).dot_product_with(incEns(0, 0));
  // check that eigen vectors are not zeros
  EXPECT(n0 > 0);
  double tol = 20*n0*DBL_EPSILON;

  // check that eig[ieig>0] are orthogonal to eig[0]
  for (int i = 1; i < neig; ++i) {
    double n = incEns(0, 0).dot_product_with(incEns(0, i));
    EXPECT(n < tol);
  }

  // try the second interface for modulateIncrement
  // this checks that modulation of a single column @geom.begin() works as well
  IncrementSet_ incEns2(Test_::resol(), Test_::ctlvars(), times, oops::mpi::myself(),
                        std::vector<int>{0});
  incEns2[0] = dx1[0];

  Eigen::MatrixXd modInc = vertloc.modulateIncrement(incEns2, geometry.begin(), 0);
  Eigen::MatrixXd modIncInner = modInc.transpose()*modInc;
  // modIncInner should be a diagonal matrix
  for (int i = 1; i < neig; ++i) {
    EXPECT(modIncInner(0, i) < modIncInner(0, 0)*DBL_EPSILON);
  }
}

// =============================================================================

template <typename MODEL>
class VerticalLocEV : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~VerticalLocEV() = default;

 private:
  std::string testid() const override {return "test::VerticalLocEV<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("generic/VerticalLocEV/testVerticalLocEV")
      { testVerticalLocEV<MODEL>(); });
  }

  void clear() const override {}
};

// =============================================================================

}  // namespace test

#endif  // TEST_GENERIC_VERTICALLOCEV_H_
