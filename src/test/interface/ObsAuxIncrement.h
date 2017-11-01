/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_OBSAUXINCREMENT_H_
#define TEST_INTERFACE_OBSAUXINCREMENT_H_

#include <iostream>
#include <string>
#include <cmath>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxCovariance.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"
#include "util/dot_product.h"

namespace test {

// =============================================================================

template <typename MODEL> class ObsAuxIncrementFixture : private boost::noncopyable {
  typedef oops::ObsAuxCovariance<MODEL> Covariance_;
  typedef oops::ObsAuxControl<MODEL>    ObsAux_;
  typedef oops::ObsAuxIncrement<MODEL>  AuxIncr_;

 public:
  static const eckit::Configuration & config()     {return *getInstance().conf_;}
  static const Covariance_  & covariance() {return *getInstance().covar_;}

 private:
  static ObsAuxIncrementFixture<MODEL>& getInstance() {
    static ObsAuxIncrementFixture<MODEL> theObsAuxIncrementFixture;
    return theObsAuxIncrementFixture;
  }

  ObsAuxIncrementFixture<MODEL>() {
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "ObsBiasCovariance"));
    covar_.reset(new Covariance_(*conf_));
  }

  ~ObsAuxIncrementFixture<MODEL>() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
  boost::scoped_ptr<const Covariance_>  covar_;
};

// =============================================================================

template <typename MODEL> void testObsAuxIncrementConstructor() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx(Test_::config());

  BOOST_CHECK_EQUAL(dx.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementCopyConstructor() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

  AuxIncr_ dx2(dx1);
  BOOST_CHECK(dx2.norm() > 0.0);
  BOOST_CHECK_EQUAL(dx2.norm(), dx1.norm());

// Check that the copy is equal to the original
  dx2 -= dx1;
  BOOST_CHECK_EQUAL(dx2.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementChangeRes() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

  AuxIncr_ dx2(dx1, Test_::config());
  BOOST_CHECK(dx2.norm() > 0.0);
  BOOST_CHECK_EQUAL(dx2.norm(), dx1.norm());

// Check that the copy is equal to the original
  dx2 -= dx1;
  BOOST_CHECK_EQUAL(dx2.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementTriangle() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx2);

// test triangle inequality
  double dot1 = dx1.norm();
  BOOST_CHECK(dot1 > 0.0);

  double dot2 = dx2.norm();
  BOOST_CHECK(dot2 > 0.0);

  dx2 += dx1;
  double dot3 = dx2.norm();
  BOOST_CHECK(dot3 > 0.0);

  BOOST_CHECK(dot3 <= dot1 + dot2);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementOpPlusEq() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(dx1);

// test *= and +=
  dx2 += dx1;
  dx1 *= 2.0;

  dx2 -= dx1;
  BOOST_CHECK_SMALL(dx2.norm(), 1e-8);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementDotProduct() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx2);

// test symmetry of dot product
  double zz1 = dot_product(dx1, dx2);
  double zz2 = dot_product(dx2, dx1);

  BOOST_CHECK_EQUAL(zz1, zz2);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementZero() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx);
  BOOST_CHECK(dx.norm() > 0.0);

// test zero
  dx->zero();
  BOOST_CHECK_EQUAL(dx.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementAxpy() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

// test axpy
  AuxIncr_ dx2(dx1);
  dx2.axpy(2.0, dx1);

  dx2 -= dx1;
  dx2 -= dx1;
  dx2 -= dx1;

  BOOST_CHECK_SMALL(dx2.norm(), 1e-8);
}

// =============================================================================

template <typename MODEL> class ObsAuxIncrement : public oops::Test {
 public:
  ObsAuxIncrement() {}
  virtual ~ObsAuxIncrement() {}
 private:
  std::string testid() const {return "test::ObsAuxIncrement<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ObsAuxIncrement");

    ts->add(BOOST_TEST_CASE(&testObsAuxIncrementConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testObsAuxIncrementCopyConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testObsAuxIncrementChangeRes<MODEL>));
    ts->add(BOOST_TEST_CASE(&testObsAuxIncrementTriangle<MODEL>));
    ts->add(BOOST_TEST_CASE(&testObsAuxIncrementOpPlusEq<MODEL>));
    ts->add(BOOST_TEST_CASE(&testObsAuxIncrementDotProduct<MODEL>));
    ts->add(BOOST_TEST_CASE(&testObsAuxIncrementAxpy<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSAUXINCREMENT_H_
