/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_INCREMENT_H_
#define TEST_INTERFACE_INCREMENT_H_

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "test/TestEnvironment.h"


namespace test {

// =============================================================================

template <typename MODEL> class IncrementFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>       Geometry_;

 public:
  static const Geometry_       & resol()   {return *getInstance().resol_;}
  static const oops::Variables & ctlvars() {return *getInstance().ctlvars_;}
  static const util::DateTime  & time()    {return *getInstance().time_;}

 private:
  static IncrementFixture<MODEL>& getInstance() {
    static IncrementFixture<MODEL> theIncrementFixture;
    return theIncrementFixture;
  }

  IncrementFixture<MODEL>() {
//  Setup a geometry
    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::comm()));

    const eckit::LocalConfiguration varConfig(TestEnvironment::config(), "Variables");
    ctlvars_.reset(new oops::Variables(varConfig));

    time_.reset(new util::DateTime(TestEnvironment::config().getString("TestDate")));
  }

  ~IncrementFixture<MODEL>() {}

  std::unique_ptr<Geometry_>       resol_;
  std::unique_ptr<oops::Variables> ctlvars_;
  std::unique_ptr<util::DateTime>  time_;
};

// =============================================================================

template <typename MODEL> void testIncrementConstructor() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());

  EXPECT(dx.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementCopyConstructor() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();
  EXPECT(dx1.norm() > 0.0);

  Increment_ dx2(dx1);
  EXPECT(dx2.norm() > 0.0);

// Check that the copy is equal to the original
  dx2 -= dx1;
  EXPECT(dx2.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementTriangle() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();
  Increment_ dx2(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx2.random();

// test triangle inequality
  double dot1 = dx1.norm();
  EXPECT(dot1 > 0.0);

  double dot2 = dx2.norm();
  EXPECT(dot2 > 0.0);

  dx2 += dx1;
  double dot3 = dx2.norm();
  EXPECT(dot3 > 0.0);

  EXPECT(dot3 <= dot1 + dot2);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementOpPlusEq() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();
  Increment_ dx2(dx1);

// test *= and +=
  dx2 += dx1;
  dx1 *= 2.0;

  dx2 -= dx1;
  EXPECT(dx2.norm()< 1e-8);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementDotProduct() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();
  Increment_ dx2(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx2.random();

// test symmetry of dot product
  double zz1 = dot_product(dx1, dx2);
  double zz2 = dot_product(dx2, dx1);

  EXPECT(zz1 == zz2);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementZero() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx.random();
  EXPECT(dx.norm() > 0.0);

// test zero
  dx->zero();
  EXPECT(dx.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementAxpy() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();

// test axpy
  Increment_ dx2(dx1);
  dx2.axpy(2.0, dx1);

  dx2 -= dx1;
  dx2 -= dx1;
  dx2 -= dx1;

  EXPECT(dx2.norm()< 1e-8);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementSerialize() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

// Create two random increments
  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();

  util::DateTime tt(Test_::time() + util::Duration("PT15H"));
  Increment_ dx2(Test_::resol(), Test_::ctlvars(), tt);

// Test serialize-deserialize
  std::vector<double> vect;
  dx1.serialize(vect);

  size_t index = 0;
  dx2.deserialize(vect, index);

  EXPECT(dx1.norm() > 0.0);
  EXPECT(dx2.norm() > 0.0);
  EXPECT(dx2.validTime() == Test_::time());

  dx2 -= dx1;
  EXPECT(dx2.norm() == 0.0);
}

// =============================================================================

template <typename MODEL>
class Increment : public oops::Test {
 public:
  Increment() {}
  virtual ~Increment() {}

 private:
  std::string testid() const {return "test::Increment<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/Increment/testIncrementConstructor")
      { testIncrementConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementCopyConstructor")
      { testIncrementCopyConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementTriangle")
      { testIncrementTriangle<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementOpPlusEq")
      { testIncrementOpPlusEq<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementDotProduct")
      { testIncrementDotProduct<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementAxpy")
      { testIncrementAxpy<MODEL>(); });
//    ts.emplace_back(CASE("interface/Increment/testIncrementSerialize")
//      { testIncrementSerialize<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_INCREMENT_H_
