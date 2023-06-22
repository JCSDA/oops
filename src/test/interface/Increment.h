/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020 UCAR.
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

#include "atlas/field.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

namespace test {

// =============================================================================

template <typename MODEL> class IncrementFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>       Geometry_;

 public:
  static const Geometry_            & resol()     {return *getInstance().resol_;}
  static const oops::Variables      & ctlvars()   {return *getInstance().ctlvars_;}
  static const util::DateTime       & time()      {return *getInstance().time_;}
  static const double               & tolerance() {return getInstance().tolerance_;}
  static const int                  & skipAtlas() {return getInstance().skipAtlas_;}
  static const eckit::Configuration & test()      {return *getInstance().test_;}
  static void reset() {
    getInstance().time_.reset();
    getInstance().ctlvars_.reset();
    getInstance().test_.reset();
    getInstance().resol_.reset();
  }

 private:
  static IncrementFixture<MODEL>& getInstance() {
    static IncrementFixture<MODEL> theIncrementFixture;
    return theIncrementFixture;
  }

  IncrementFixture<MODEL>() {
//  Setup a geometry
    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::world()));

    ctlvars_.reset(new oops::Variables(TestEnvironment::config(), "inc variables"));

    const double tol_default = 1e-8;
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "increment test"));
    time_.reset(new util::DateTime(test_->getString("date")));
    tolerance_ = test_->getDouble("tolerance", tol_default);
    if (tolerance_ > tol_default) {
      oops::Log::warning() <<
        "Warning: Increment norm tolerance greater than 1e-8 "
        "may not be suitable for certain solvers." <<
        std::endl; }
    skipAtlas_ = test_->getInt("skip atlas", 0);
  }

  ~IncrementFixture<MODEL>() {}

  std::unique_ptr<Geometry_>       resol_;
  std::unique_ptr<oops::Variables> ctlvars_;
  std::unique_ptr<const eckit::LocalConfiguration> test_;
  double                           tolerance_;
  int                              skipAtlas_;
  std::unique_ptr<util::DateTime>  time_;
  std::unique_ptr<bool> skipAccumTest_;
  std::unique_ptr<bool> skipDiffTest_;
  std::unique_ptr<bool> skipRmsByLevelTest_;
};

// =============================================================================
/// \brief tests Increment constructor and print method
template <typename MODEL> void testIncrementConstructor() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());
  oops::Log::test() << "Printing zero increment: " << dx << std::endl;

  EXPECT(dx.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementCopyConstructor() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();
  oops::Log::test() << "Printing random increment: " << dx1 << std::endl;

  EXPECT(dx1.norm() > 0.0);

  Increment_ dx2(dx1);
  EXPECT(dx2.norm() > 0.0);

// Check that the copy is equal to the original
  dx2 -= dx1;
  EXPECT(dx2.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementCopyBoolConstructor() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();
  EXPECT(dx1.norm() > 0.0);

  // Test with copy set to true
  Increment_ dx2(dx1, true);
  EXPECT(dx2.norm() == dx1.norm());

  // Test with copy set to false
  Increment_ dx3(dx1, false);
  EXPECT(dx3.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementChangeResConstructor() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  // For now this is just a copy and not a change res, would require config changes
  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  Increment_ dx2(Test_::resol(), dx1);

  // Check they are same. Should be replaced with change res and check they are different
  EXPECT(dx2.norm() == dx1.norm());
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
  EXPECT(dx2.norm() < Test_::tolerance());
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
  dx.zero();
  EXPECT(dx.norm() == 0.0);
  EXPECT(dx.validTime() == Test_::time());

// Create a new time one hour in the future
  util::Duration onehour(3600);
  util::DateTime newTime(Test_::time().toString());
  newTime+=onehour;

// Confirm zero with setting new datetime works
  dx.random();
  dx.zero(newTime);
  EXPECT(dx.norm() == 0.0);
  EXPECT(dx.validTime() == newTime);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementDirac() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

// Option to skip test
  bool skipTest = Test_::test().getBool("skip dirac test", false);
  if (skipTest) {
    oops::Log::warning() << "Skipping Increment.dirac test";
    return;
  }

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());

// test dirac
  dx.dirac(Test_::test().getSubConfiguration("dirac"));
  Increment_ dx2_minus_dx(dx);
  dx2_minus_dx.schur_product_with(dx);
  dx2_minus_dx -= dx;
  EXPECT(dx.norm() > 0.0);
  EXPECT(dx2_minus_dx.norm() == 0.0);
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

  EXPECT(dx2.norm() < Test_::tolerance());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementAccum() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;
  typedef oops::State<MODEL>        State_;

  // Option to skip test
  bool skipTest = Test_::test().getBool("skip accum test", false);
  if (skipTest) {
    oops::Log::warning() << "Skipping Increment.accum test";
    return;
  }

  // Create two different random increments
  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();

  Increment_ dx2(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx2.random();

  Increment_ diff(dx1);
  diff -= dx2;
  EXPECT(diff.norm() != 0.0);

  // Create a state that is equal to dx2
  State_ x(Test_::resol(), Test_::ctlvars(), Test_::time());
  x.zero();
  x += dx2;

  // Create copy of dx1 to test against axpy
  Increment_ dx3(dx1);

  // Test accum
  dx1.accumul(2.0, x);

  // Use axpy for reference
  dx3.axpy(2.0, dx2);

  // Check axpy did something
  diff = dx3;
  diff -= dx2;
  EXPECT(diff.norm() != 0.0);

  // Check accumul matches axpy
  diff = dx1;
  diff -= dx3;
  EXPECT(diff.norm() == 0.0);
}

// -----------------------------------------------------------------------------
template <typename MODEL> void testRmsByLevel() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  // Option to skip test
  const bool skipTest = Test_::test().getBool("skip rms by level test", false);
  if (skipTest) {
    oops::Log::warning() << "Skipping Increment.rmsByLevel test";
    return;
  } else {
      Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
      dx1.ones();
      std::vector<double> vec = dx1.rmsByLevel(dx1.variables()[0]);
      std::vector<double> referenceVec(vec.size(), 1.0);
      EXPECT(vec == referenceVec);
  }
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
  EXPECT(vect.size() == dx1.serialSize());

  size_t index = 0;
  dx2.deserialize(vect, index);
  EXPECT(index == dx1.serialSize());
  EXPECT(index == dx2.serialSize());

  dx1.serialize(vect);
  EXPECT(vect.size() == dx1.serialSize() * 2);

  if (dx1.serialSize() > 0) {  // until all models have implemented serialize
    EXPECT(dx1.norm() > 0.0);
    EXPECT(dx2.norm() > 0.0);
    EXPECT(dx2.validTime() == Test_::time());

    dx2 -= dx1;
    EXPECT(dx2.norm() == 0.0);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementDiff() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;
  typedef oops::State<MODEL>        State_;

  // Option to skip test
  bool skipTest = Test_::test().getBool("skip diff test", false);
  if (skipTest) {
    oops::Log::warning() << "Skipping Increment.diff test";
    return;
  }

  // Create two states to diff
  State_ x1(Test_::resol(), Test_::ctlvars(), Test_::time());
  State_ x2(Test_::resol(), Test_::ctlvars(), Test_::time());

  // Create two different random increments
  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();

  Increment_ dx2(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx2.random();

  // Create some different increments
  Increment_ diff1(dx1);
  Increment_ diff2(Test_::resol(), Test_::ctlvars(), Test_::time());
  diff1 -= dx2;
  EXPECT(diff1.norm() != 0.0);

  // Fill states with some different random values
  x1.zero();
  x1+=dx1;
  x2.zero();
  x2+=dx2;

  // Use difference of states to compute difference
  diff2.diff(x1, x2);

  // Compare difference with -= operator
  EXPECT(diff1.norm() == diff2.norm());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementTime() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  // Create an increment
  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());

  // Move increment time forward by one hour
  util::Duration onehour(3600);
  dx.updateTime(onehour);

  // Confirm new valid time and that validTime function works
  util::DateTime newTime(Test_::time().toString());
  newTime+=onehour;
  EXPECT(dx.validTime() == newTime);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementSchur() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  // Create two identical random increments
  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();
  Increment_ dx2(dx1);
  EXPECT(dx1.norm() == dx2.norm());

  // Compute schur product
  dx1.schur_product_with(dx2);

  // Checks
  EXPECT(dx1.norm() != dx2.norm());

  // Pass through zeros
  dx1.random();
  dx2.zero();
  dx1.schur_product_with(dx2);
  EXPECT(dx1.norm() == 0.0);

  // Pass through ones
  dx1.random();
  Increment_ dx1in(dx1);
  dx2.ones();
  dx1.schur_product_with(dx2);
  EXPECT(dx1.norm() == dx1in.norm());
}

// =============================================================================

template <typename MODEL>
class Increment : public oops::Test {
 public:
  Increment() {}
  virtual ~Increment() {IncrementFixture<MODEL>::reset();}

 private:
  std::string testid() const override {return "test::Increment<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/Increment/testIncrementConstructor")
      { testIncrementConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementCopyConstructor")
      { testIncrementCopyConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementCopyBoolConstructor")
      { testIncrementCopyBoolConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementChangeResConstructor")
      { testIncrementChangeResConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/rmsByLevel")
      { testRmsByLevel<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementTriangle")
      { testIncrementTriangle<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementOpPlusEq")
      { testIncrementOpPlusEq<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementDotProduct")
      { testIncrementDotProduct<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementAxpy")
      { testIncrementAxpy<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementAccum")
      { testIncrementAccum<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementDiff")
      { testIncrementDiff<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementZero")
      { testIncrementZero<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementDirac")
      { testIncrementDirac<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementTime")
      { testIncrementTime<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementSchur")
      { testIncrementSchur<MODEL>(); });
    ts.emplace_back(CASE("interface/Increment/testIncrementSerialize")
      { testIncrementSerialize<MODEL>(); });
  }

  void clear() const override {}
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_INCREMENT_H_
