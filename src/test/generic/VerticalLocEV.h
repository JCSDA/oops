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
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementEnsemble4D.h"
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
  typedef oops::IncrementEnsemble4D<MODEL>        IncrementEnsemble_;
  typedef oops::State<MODEL>              State_;

  const Geometry_ & geometry = Test_::resol();
  eckit::LocalConfiguration vertlocconf(TestEnvironment::config(), "vertical localization");

  // make an empty state vector to be used to intitialize VerticalLocEV_
  State_ x(Test_::resol(), Test_::ctlvars(), Test_::time());

  VerticalLocEV_ vertloc(vertlocconf, x, x.variables());
  oops::Log::test() << "Number of eigenvalues used in VerticalLoc: " << vertloc.neig() << std::endl;

  //--- check for expected number of eigen modes
  int nEigExpected = TestEnvironment::config().getInt("expected neig");
  int neig = vertloc.neig();
  oops::Log::debug() << "Expected number of eigen modes: " << nEigExpected << std::endl;
  oops::Log::debug() << "Actual number of eigen modes: " << neig << std::endl;
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
  IncrementEnsemble_ incEns(Test_::resol(), Test_::ctlvars(), times, neig);
  // set incEns to zero
  for (int i = 1; i < neig; ++i) {incEns[i][0].zero();}
  for (int i = 1; i < neig; ++i) {
    double n = incEns[0][0].dot_product_with(incEns[i][0]);
    EXPECT(n < 100*DBL_EPSILON);
  }

  // make increment of ones and check that it is indeed full of ones
  dx1.ones();
  dx2.random();
  double normRandom = dx2.dot_product_with(dx2);
  dx2.schur_product_with(dx1);
  EXPECT(std::abs(dx2.dot_product_with(dx2)-normRandom) < normRandom*DBL_EPSILON);
  oops::Log::debug() << "Increment ones()" << dx1 << std::endl;

  // modulate increments
  vertloc.modulateIncrement(dx1, incEns);

  // check the orthogonality condition
  double n0 = incEns[0][0].dot_product_with(incEns[0][0]);
  oops::Log::debug() << "dot product 0: " <<  n0 << std::endl;
  // check that eigen vectors are not zeros
  EXPECT(n0 > 0);
  double tol = 20*n0*DBL_EPSILON;
  oops::Log::debug() << "tolerance :" << tol << std::endl;

  // check that eig[ieig>0] are orthogonal to eig[0]
  for (int i = 1; i < neig; ++i) {
    double n = incEns[0][0].dot_product_with(incEns[i][0]);
    oops::Log::debug() << "dot product " << i << ": " <<  n << std::endl;
    EXPECT(n < tol);
  }

  // try the second interface for modulateIncrement
  // this checks that modulation of a single column @geom.begin() works as well
  IncrementEnsemble_ incEns2(Test_::resol(), Test_::ctlvars(), times, 1);
  incEns2[0] = dx1;

  Eigen::MatrixXd modInc = vertloc.modulateIncrement(incEns2, geometry.begin(), 0);
  Eigen::MatrixXd modIncInner = modInc.transpose()*modInc;
  oops::Log::debug() << "modInc'*modInc" << modIncInner << std::endl;
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
