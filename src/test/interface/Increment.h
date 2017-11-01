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

#include <iostream>
#include <string>
#include <cmath>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/runs/Test.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "test/TestEnvironment.h"
#include "util/DateTime.h"
#include "util/dot_product.h"

namespace test {

// =============================================================================

template <typename MODEL> class IncrementFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>       Geometry_;
  typedef oops::State<MODEL>          State_;
  typedef oops::Variables<MODEL>      Variables_;

 public:
  static const Geometry_      & resol()      {return *getInstance().resol_;}
  static const Variables_     & ctlvars()    {return *getInstance().ctlvars_;}
  static const util::DateTime & time()       {return *getInstance().time_;}

 private:
  static IncrementFixture<MODEL>& getInstance() {
    static IncrementFixture<MODEL> theIncrementFixture;
    return theIncrementFixture;
  }

  IncrementFixture<MODEL>() {
//  Setup a geometry
    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));

    const eckit::LocalConfiguration varConfig(TestEnvironment::config(), "Variables");
    ctlvars_.reset(new Variables_(varConfig));

//  Setup reference state
    const eckit::LocalConfiguration fgconf(TestEnvironment::config(), "State");
    State_ xx(*resol_, fgconf);

    time_.reset(new util::DateTime(xx.validTime()));
  }

  ~IncrementFixture<MODEL>() {}

  boost::scoped_ptr<Geometry_>      resol_;
  boost::scoped_ptr<Variables_>     ctlvars_;
  boost::scoped_ptr<util::DateTime> time_;
};

// =============================================================================

template <typename MODEL> void testIncrementConstructor() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());

  BOOST_CHECK_EQUAL(dx.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementCopyConstructor() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();
  BOOST_CHECK(dx1.norm() > 0.0);

  Increment_ dx2(dx1);
  BOOST_CHECK(dx2.norm() > 0.0);

// Check that the copy is equal to the original
  dx2 -= dx1;
  BOOST_CHECK_EQUAL(dx2.norm(), 0.0);
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
  BOOST_CHECK(dot1 > 0.0);

  double dot2 = dx2.norm();
  BOOST_CHECK(dot2 > 0.0);

  dx2 += dx1;
  double dot3 = dx2.norm();
  BOOST_CHECK(dot3 > 0.0);

  BOOST_CHECK(dot3 < dot1 + dot2);
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
  BOOST_CHECK_SMALL(dx2.norm(), 1e-8);
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

  BOOST_CHECK_EQUAL(zz1, zz2);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testIncrementZero() {
  typedef IncrementFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx.random();
  BOOST_CHECK(dx.norm() > 0.0);

// test zero
  dx->zero();
  BOOST_CHECK_EQUAL(dx.norm(), 0.0);
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

  BOOST_CHECK_SMALL(dx2.norm(), 1e-8);
}

// =============================================================================

template <typename MODEL> class Increment : public oops::Test {
 public:
  Increment() {}
  virtual ~Increment() {}
 private:
  std::string testid() const {return "test::Increment<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/Increment");

    ts->add(BOOST_TEST_CASE(&testIncrementConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testIncrementCopyConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testIncrementTriangle<MODEL>));
    ts->add(BOOST_TEST_CASE(&testIncrementOpPlusEq<MODEL>));
    ts->add(BOOST_TEST_CASE(&testIncrementDotProduct<MODEL>));
    ts->add(BOOST_TEST_CASE(&testIncrementAxpy<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_INCREMENT_H_
