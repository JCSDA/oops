/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_MODELAUXINCREMENT_H_
#define TEST_INTERFACE_MODELAUXINCREMENT_H_

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
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxCovariance.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"
#include "util/dot_product.h"

namespace test {

// =============================================================================

template <typename MODEL> class ModelAuxIncrementFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>           Geometry_;
  typedef oops::ModelAuxCovariance<MODEL> Covariance_;
  typedef oops::ModelAuxControl<MODEL>    ModelAux_;
  typedef oops::ModelAuxIncrement<MODEL>  AuxIncr_;

 public:
  static const eckit::Configuration & config()     {return *getInstance().conf_;}
  static const Covariance_  & covariance() {return *getInstance().covar_;}
  static const Geometry_    & resol()      {return *getInstance().resol_;}

 private:
  static ModelAuxIncrementFixture<MODEL>& getInstance() {
    static ModelAuxIncrementFixture<MODEL> theModelAuxIncrementFixture;
    return theModelAuxIncrementFixture;
  }

  ModelAuxIncrementFixture<MODEL>() {
//  Setup a geometry
    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));

//  Setup a covariance matrix
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "ModelBiasCovariance"));
    covar_.reset(new Covariance_(*conf_, *resol_));
  }

  ~ModelAuxIncrementFixture<MODEL>() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
  boost::scoped_ptr<const Geometry_>    resol_;
  boost::scoped_ptr<const Covariance_>  covar_;
};

// =============================================================================

template <typename MODEL> void testModelAuxIncrementConstructor() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx(Test_::resol(), Test_::config());

  BOOST_CHECK_EQUAL(dx.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementCopyConstructor() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

  AuxIncr_ dx2(dx1);
  BOOST_CHECK(dx2.norm() > 0.0);
  BOOST_CHECK_EQUAL(dx2.norm(), dx1.norm());

// Check that the copy is equal to the original
  dx2 -= dx1;
  BOOST_CHECK_EQUAL(dx2.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementChangeRes() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

  AuxIncr_ dx2(dx1, Test_::config());
  BOOST_CHECK(dx2.norm() > 0.0);
  BOOST_CHECK_EQUAL(dx2.norm(), dx1.norm());

// Check that the copy is equal to the original
  dx2 -= dx1;
  BOOST_CHECK_EQUAL(dx2.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementTriangle() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx2);

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

template <typename MODEL> void testModelAuxIncrementOpPlusEq() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(dx1);

// test *= and +=
  dx2 += dx1;
  dx1 *= 2.0;

  dx2 -= dx1;
  BOOST_CHECK_SMALL(dx2.norm(), 1e-8);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementDotProduct() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx2);

// test symmetry of dot product
  double zz1 = dot_product(dx1, dx2);
  double zz2 = dot_product(dx2, dx1);

  BOOST_CHECK_EQUAL(zz1, zz2);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementZero() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx);
  BOOST_CHECK(dx.norm() > 0.0);

// test zero
  dx->zero();
  BOOST_CHECK_EQUAL(dx.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementAxpy() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

// test axpy
  AuxIncr_ dx2(dx1);
  dx2.axpy(2.0, dx1);

  dx2 -= dx1;
  dx2 -= dx1;
  dx2 -= dx1;

  BOOST_CHECK_SMALL(dx2.norm(), 1e-8);
}

// =============================================================================

template <typename MODEL> class ModelAuxIncrement : public oops::Test {
 public:
  ModelAuxIncrement() {}
  virtual ~ModelAuxIncrement() {}
 private:
  std::string testid() const {return "test::ModelAuxIncrement<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/ModelAuxIncrement");

    ts->add(BOOST_TEST_CASE(&testModelAuxIncrementConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testModelAuxIncrementCopyConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testModelAuxIncrementChangeRes<MODEL>));
    ts->add(BOOST_TEST_CASE(&testModelAuxIncrementTriangle<MODEL>));
    ts->add(BOOST_TEST_CASE(&testModelAuxIncrementOpPlusEq<MODEL>));
    ts->add(BOOST_TEST_CASE(&testModelAuxIncrementDotProduct<MODEL>));
    ts->add(BOOST_TEST_CASE(&testModelAuxIncrementAxpy<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_MODELAUXINCREMENT_H_
