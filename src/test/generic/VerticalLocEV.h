/*
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_GENERIC_VERTICALLOCEV_H_
#define TEST_GENERIC_VERTICALLOCEV_H_

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/base/Variables.h"
#include "oops/generic/VerticalLocEV.h"
#include "oops/interface/Geometry.h"
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
  typedef oops::Increment4D<MODEL>        Increment4D_;
  typedef oops::IncrementEnsemble<MODEL>  IncrementEnsemble_;
  typedef oops::VerticalLocEV<MODEL>      VerticalLocEV_;

  const Geometry_ & geometry = Test_::resol();
  eckit::LocalConfiguration vertlocconf(TestEnvironment::config(), "vertical localization");
  VerticalLocEV_ vertloc(geometry, vertlocconf);
  oops::Log::test() << "Number of eigenvalues used in VerticalLoc: " << vertloc.neig() << std::endl;

  // check for expected number of eigen modes
  int nEigExpected = TestEnvironment::config().getDouble("expected neig");
  int neig = vertloc.neig();
  oops::Log::debug() << "Expected number of eigen modes: " << nEigExpected << std::endl;
  oops::Log::debug() << "Actual number of eigen modes: " << neig << std::endl;
  EXPECT(nEigExpected == neig);

  // check that truncation and rescaling was done correctley
  EXPECT(vertloc.testTruncateEvecs());

/* TODO(Issue #828) finish unit test for modulated product
  // check modulation
  Increment4D_ dx(Test_::resol(), Test_::ctlvars(), times);
  dx.ones();

  oops::Log::test() << "Increment ones()" << dx << std::endl;
  oops::Log::test() << "norm(dx):" << dx.dot_product_with(dx) << std::endl;

  IncrementEnsemble_ incEns(Test_::resol(), Test_::ctlvars(), times, neig);
  vertloc.modulateIncrement(dx, incEns);

  oops::Log::debug() << "dot product 0: " <<  incEns[0][0].dot_product_with(incEns[0][0]) << std::endl;
  oops::Log::debug() << "dot product 1: " <<  incEns[1][0].dot_product_with(incEns[1][0]) << std::endl;

//  Eigen::MatrixXd modulateIncrement(const IncrementEnsemble_ &,
//                                    const GeometryIterator_ &, size_t) const;
*/
}

// =============================================================================

template <typename MODEL>
class VerticalLocEV : public oops::Test {
 public:
  VerticalLocEV() = default;
  virtual ~VerticalLocEV() = default;

 private:
  std::string testid() const {return "test::VerticalLocEV<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("generic/VerticalLocEV/testVerticalLocEV")
      { testVerticalLocEV<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_GENERIC_VERTICALLOCEV_H_
