/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_LINEARVARIABLECHANGE_H_
#define TEST_INTERFACE_LINEARVARIABLECHANGE_H_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/generic/instantiateLinearVariableChangeFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> class LinearVariableChangeFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>                  Geometry_;
  typedef oops::State<MODEL>                     State_;
  typedef util::DateTime                         DateTime_;

 public:
  static std::vector<eckit::LocalConfiguration>
                                     & linvarchgconfs()      {return getInstance().linvarchgconfs_;}
  static const State_                & xx()                  {return *getInstance().xx_;}
  static const Geometry_             & resol()               {return *getInstance().resol_;}
  static const DateTime_             & time()                {return *getInstance().time_;}

 private:
  static LinearVariableChangeFixture<MODEL>& getInstance() {
    static LinearVariableChangeFixture<MODEL> theLinearVariableChangeFixture;
    return theLinearVariableChangeFixture;
  }

  LinearVariableChangeFixture<MODEL>() {
    oops::instantiateLinearVariableChangeFactory<MODEL>();

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));

    const eckit::LocalConfiguration fgconf(TestEnvironment::config(), "State");
    xx_.reset(new State_(*resol_, fgconf));

    time_.reset(new util::DateTime(xx_->validTime()));

    TestEnvironment::config().get("LinearVariableChangeTests", linvarchgconfs_);
  }

  ~LinearVariableChangeFixture<MODEL>() {}

  std::vector<eckit::LocalConfiguration>             linvarchgconfs_;
  boost::scoped_ptr<const State_ >                   xx_;
  boost::scoped_ptr<const Geometry_>                 resol_;
  boost::scoped_ptr<const util::DateTime>            time_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVariableChangeZero() {
  typedef LinearVariableChangeFixture<MODEL>       Test_;
  typedef oops::Increment<MODEL>                   Increment_;
  typedef oops::LinearVariableChangeBase<MODEL>    LinearVariableChange_;
  typedef oops::LinearVariableChangeFactory<MODEL> LinearVariableChangeFactory_;

  for (std::size_t jj = 0; jj < Test_::linvarchgconfs().size(); ++jj) {
    eckit::LocalConfiguration varinconf(Test_::linvarchgconfs()[jj], "inputVariables");
    eckit::LocalConfiguration varoutconf(Test_::linvarchgconfs()[jj], "outputVariables");
    oops::Variables varin(varinconf);
    oops::Variables varout(varoutconf);

    boost::scoped_ptr<LinearVariableChange_> changevar(LinearVariableChangeFactory_::create(
                                      Test_::xx(), Test_::xx(),
                                      Test_::resol(), Test_::linvarchgconfs()[jj]));

    Increment_   dxin(Test_::resol(), varin,  Test_::time());
    Increment_ KTdxin(Test_::resol(), varout, Test_::time());
    Increment_  dxout(Test_::resol(), varout, Test_::time());
    Increment_ Kdxout(Test_::resol(), varin,  Test_::time());

    // dxout = 0, check if K.dxout = 0
    dxout.zero();
    changevar->multiply(dxout, Kdxout);
    BOOST_CHECK_EQUAL(Kdxout.norm(), 0.0);

    // dxin = 0, check if K^T.dxin = 0
    dxin.zero();
    changevar->multiplyAD(dxin, KTdxin);
    BOOST_CHECK_EQUAL(KTdxin.norm(), 0.0);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVariableChangeInverse() {
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVariableChangeAdjoint() {
  typedef LinearVariableChangeFixture<MODEL>       Test_;
  typedef oops::Increment<MODEL>                   Increment_;
  typedef oops::LinearVariableChangeBase<MODEL>    LinearVariableChange_;
  typedef oops::LinearVariableChangeFactory<MODEL> LinearVariableChangeFactory_;

  for (std::size_t jj = 0; jj < Test_::linvarchgconfs().size(); ++jj) {
    eckit::LocalConfiguration varinconf(Test_::linvarchgconfs()[jj], "inputVariables");
    eckit::LocalConfiguration varoutconf(Test_::linvarchgconfs()[jj], "outputVariables");
    oops::Variables varin(varinconf);
    oops::Variables varout(varoutconf);

    boost::scoped_ptr<LinearVariableChange_> changevar(LinearVariableChangeFactory_::create(
                                      Test_::xx(), Test_::xx(),
                                      Test_::resol(), Test_::linvarchgconfs()[jj]));

    Increment_   dxin(Test_::resol(), varin,  Test_::time());
    Increment_ KTdxin(Test_::resol(), varout, Test_::time());
    Increment_  dxout(Test_::resol(), varout, Test_::time());
    Increment_ Kdxout(Test_::resol(), varin,  Test_::time());

    dxin.random();
    dxout.random();

    Increment_  dxin0(dxin);
    Increment_  dxout0(dxout);

    changevar->multiply(dxout, Kdxout);
    changevar->multiplyAD(dxin, KTdxin);

    // zz1 = <Kdxout,dxin>
    const double zz1 = dot_product(Kdxout, dxin0);
    // zz2 = <dxout,KTdxin>
    const double zz2 = dot_product(dxout0, KTdxin);

    // Print
    oops::Log::info() << "<dxout,KTdxin>-<Kdxout,dxin>/<dxout,KTdxin>="
                      << (zz1-zz2)/zz1 << std::endl;
    oops::Log::info() << "<dxout,KTdxin>-<Kdxout,dxin>/<Kdxout,dxin>="
                      << (zz1-zz2)/zz2 << std::endl;
    const double tol = 1e-8;
    BOOST_CHECK_CLOSE(zz1, zz2, tol);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> class LinearVariableChange : public oops::Test {
 public:
  LinearVariableChange() {}
  virtual ~LinearVariableChange() {}
 private:
  std::string testid() const {return "test::LinearVariableChange<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/LinearVariableChange");

    ts->add(BOOST_TEST_CASE(&testLinearVariableChangeZero<MODEL>));
    ts->add(BOOST_TEST_CASE(&testLinearVariableChangeAdjoint<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LINEARVARIABLECHANGE_H_
