/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_VARIABLECHANGE_H_
#define TEST_INTERFACE_VARIABLECHANGE_H_

#include <cmath>
#include <iostream>
#include <string>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/VariableChangeBase.h"
#include "oops/generic/instantiateVariableChangeFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

namespace test {

// =============================================================================

template <typename MODEL> class VariableChangeFixture : private boost::noncopyable {
  typedef oops::VariableChangeBase<MODEL>  VariableChange_;
  typedef oops::Geometry<MODEL>            Geometry_;

 public:
  static const eckit::Configuration & test()      {return *getInstance().test_;}
  static const Geometry_            & resol()     {return *getInstance().resol_;}
  static const oops::Variables      & varin()     {return *getInstance().varin_;}
  static const oops::Variables      & varout()    {return *getInstance().varout_;}
  static const util::DateTime       & time()      {return *getInstance().time_;}
  static const VariableChange_      & changevar() {return *getInstance().K_;}

 private:
  static VariableChangeFixture<MODEL>& getInstance() {
    static VariableChangeFixture<MODEL> theVariableChangeFixture;
    return theVariableChangeFixture;
  }

  VariableChangeFixture<MODEL>() {
    oops::instantiateVariableChangeFactory<MODEL>();

    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "VariableChangeTest"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));

    const eckit::LocalConfiguration fgconf(TestEnvironment::config(), "State");
    oops::State<MODEL> xx(*resol_, fgconf);

    time_.reset(new util::DateTime(xx.validTime()));

//  Setup the change of variable
    const eckit::LocalConfiguration changevarconf(TestEnvironment::config(), "VariableChangeTest");
    const eckit::LocalConfiguration varinconf(TestEnvironment::config(),
                                          "VariableChangeTest.inputVariables");
    varin_.reset(new oops::Variables(varinconf));

    const eckit::LocalConfiguration varoutconf(TestEnvironment::config(),
                                          "VariableChangeTest.outputVariables");
    varout_.reset(new oops::Variables(varoutconf));

    K_.reset(oops::VariableChangeFactory<MODEL>::create(changevarconf));
    K_->linearize(xx, *resol_);
  }

  ~VariableChangeFixture<MODEL>() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> test_;
  boost::scoped_ptr<const Geometry_>                 resol_;
  boost::scoped_ptr<const oops::Variables>           varin_;
  boost::scoped_ptr<const oops::Variables>           varout_;
  boost::scoped_ptr<const util::DateTime>            time_;
  boost::scoped_ptr<VariableChange_>                 K_;
};

// =============================================================================

template <typename MODEL> void testVariableChangeZero() {
  typedef VariableChangeFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dxin(Test_::resol(), Test_::varin(), Test_::time());
  Increment_ dxout(Test_::resol(), Test_::varout(), Test_::time());

  // dxin = 0, check if K.dxin = 0
  dxin.zero();
  dxout.zero();
  dxout = Test_::changevar().multiply(dxin);
  BOOST_CHECK_EQUAL(dxout.norm(), 0.0);

  // dxout = 0
  // test K^T.dxout = 0
  dxin.zero();
  dxout.zero();
  dxin = Test_::changevar().multiplyAD(dxout);
  BOOST_CHECK_EQUAL(dxin.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testVariableChangeInverse() {
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testVariableChangeAdjoint() {
  typedef VariableChangeFixture<MODEL>   Test_;
  typedef oops::Increment<MODEL>    Increment_;

  Increment_ dxin(Test_::resol(), Test_::varin(), Test_::time());
  Increment_ Kdxin(Test_::resol(), Test_::varout(), Test_::time());
  Increment_ dxout(Test_::resol(), Test_::varout(), Test_::time());
  Increment_ KTdxout(Test_::resol(), Test_::varin(), Test_::time());

  dxin.random();
  dxout.random();

  const double zz1 = dot_product(dxout, Test_::changevar().multiply(dxin));
  const double zz2 = dot_product(Test_::changevar().multiply(dxout), dxin);
  oops::Log::info() << "<dxout,Kdxin>-<KTdxout,dxin>/<dxout,Kdxin>="
                    <<  (zz1-zz2)/zz1 << std::endl;
  oops::Log::info() << "<dxout,Kdxin>-<KTdxout,dxin>/<KTdxout,dxin>="
                    <<  (zz1-zz2)/zz2 << std::endl;
  const double tol = 1e-8;
  BOOST_CHECK_CLOSE(zz1, zz2, tol);
}

// =============================================================================

template <typename MODEL> class VariableChange : public oops::Test {
 public:
  VariableChange() {}
  virtual ~VariableChange() {}
 private:
  std::string testid() const {return "test::VariableChange<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/VariableChange");

    ts->add(BOOST_TEST_CASE(&testVariableChangeZero<MODEL>));
    ts->add(BOOST_TEST_CASE(&testVariableChangeAdjoint<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_VARIABLECHANGE_H_
