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
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/testing/Test.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/Variables.h"
#include "oops/generic/instantiateVariableChangeFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
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
  static std::vector<eckit::LocalConfiguration> & confs() {return getInstance().confs_;}
  static const State_      & xx()               {return *getInstance().xx_;}
  static const Geometry_   & resol()            {return *getInstance().resol_;}
  static const DateTime_   & time()             {return *getInstance().time_;}

 private:
  static LinearVariableChangeFixture<MODEL>& getInstance() {
    static LinearVariableChangeFixture<MODEL> theLinearVariableChangeFixture;
    return theLinearVariableChangeFixture;
  }

  LinearVariableChangeFixture<MODEL>() {
    oops::instantiateVariableChangeFactory<MODEL>();

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::comm()));

    const eckit::LocalConfiguration fgconf(TestEnvironment::config(), "background");
    xx_.reset(new State_(*resol_, fgconf));

    time_.reset(new util::DateTime(xx_->validTime()));

    TestEnvironment::config().get("linear variable change tests", confs_);
  }

  ~LinearVariableChangeFixture<MODEL>() {}

  std::vector<eckit::LocalConfiguration>             confs_;
  std::unique_ptr<const State_ >                   xx_;
  std::unique_ptr<const Geometry_>                 resol_;
  std::unique_ptr<const util::DateTime>            time_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVariableChangeZero() {
  typedef LinearVariableChangeFixture<MODEL>       Test_;
  typedef oops::Increment<MODEL>                   Increment_;
  typedef oops::LinearVariableChangeBase<MODEL>    LinearVariableChange_;
  typedef oops::LinearVariableChangeFactory<MODEL> LinearVariableChangeFactory_;

  for (std::size_t jj = 0; jj < Test_::confs().size(); ++jj) {
    oops::Variables varin(Test_::confs()[jj], "input variables");
    oops::Variables varout(Test_::confs()[jj], "output variables");

    std::unique_ptr<LinearVariableChange_> changevar(LinearVariableChangeFactory_::create(
                                      Test_::xx(), Test_::xx(),
                                      Test_::resol(), Test_::confs()[jj]));

    Increment_   dxin(Test_::resol(), varin,  Test_::time());
    Increment_ KTdxin(Test_::resol(), varout, Test_::time());
    Increment_  dxout(Test_::resol(), varout, Test_::time());
    Increment_ Kdxout(Test_::resol(), varin,  Test_::time());

    // dxout = 0, check if K.dxout = 0
    dxout.zero();
    changevar->multiply(dxout, Kdxout);
    EXPECT(Kdxout.norm() == 0.0);

    // dxin = 0, check if K^T.dxin = 0
    dxin.zero();
    changevar->multiplyAD(dxin, KTdxin);
    EXPECT(KTdxin.norm() == 0.0);

    const bool testinverse = Test_::confs()[jj].getBool("test inverse", true);
    if (testinverse)
      {
        Increment_   KIdxin(Test_::resol(), varout, Test_::time());
        Increment_ KTIdxout(Test_::resol(), varin,  Test_::time());

        oops::Log::info() << "Doing zero test for inverse" << std::endl;
        dxout.zero();
        changevar->multiplyInverseAD(dxout, KTIdxout);
        EXPECT(KTIdxout.norm() == 0.0);

        dxin.zero();
        changevar->multiplyInverse(dxin, KIdxin);
        EXPECT(KIdxin.norm() == 0.0);
      } else {
      oops::Log::info() << "Not doing zero test for inverse" << std::endl;
    }
  }
}
// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVariableChangeAdjoint() {
  typedef LinearVariableChangeFixture<MODEL>       Test_;
  typedef oops::Increment<MODEL>                   Increment_;
  typedef oops::LinearVariableChangeBase<MODEL>    LinearVariableChange_;
  typedef oops::LinearVariableChangeFactory<MODEL> LinearVariableChangeFactory_;

  for (std::size_t jj = 0; jj < Test_::confs().size(); ++jj) {
    oops::Variables varin(Test_::confs()[jj], "input variables");
    oops::Variables varout(Test_::confs()[jj], "output variables");

    std::unique_ptr<LinearVariableChange_> changevar(LinearVariableChangeFactory_::create(
                                      Test_::xx(), Test_::xx(),
                                      Test_::resol(), Test_::confs()[jj]));

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
    double zz1 = dot_product(Kdxout, dxin0);
    // zz2 = <dxout,KTdxin>
    double zz2 = dot_product(dxout0, KTdxin);

    oops::Log::info() << "<dxout,KTdxin>-<Kdxout,dxin>/<dxout,KTdxin>="
                      << (zz1-zz2)/zz1 << std::endl;
    oops::Log::info() << "<dxout,KTdxin>-<Kdxout,dxin>/<Kdxout,dxin>="
                      << (zz1-zz2)/zz2 << std::endl;
    const double tol = 1e-10;
    EXPECT(oops::is_close(zz1, zz2, tol));
    const bool testinverse = Test_::confs()[jj].getBool("test inverse", true);
    if (testinverse)
      {
        Increment_   invKdxin(Test_::resol(), varout, Test_::time());
        Increment_ KTIdxout(Test_::resol(), varin,  Test_::time());
        oops::Log::info() << "Doing adjoint test for inverse" << std::endl;
        dxin.random();
        dxout.random();
        dxin0 = dxin;
        dxout0 = dxout;
        changevar->multiplyInverseAD(dxout, KTIdxout);
        changevar->multiplyInverse(dxin, invKdxin);
        zz1 = dot_product(KTIdxout, dxin0);
        zz2 = dot_product(dxout0, invKdxin);
        oops::Log::info() << "<dxout,KinvTdxin>-<Kinvdxout,dxin>/<dxout,KinvTdxin>="
                      << (zz1-zz2)/zz1 << std::endl;
        oops::Log::info() << "<dxout,KinvTdxin>-<Kinvdxout,dxin>/<Kinvdxout,dxin>="
                      << (zz1-zz2)/zz2 << std::endl;
        EXPECT(oops::is_close(zz1, zz2, tol));
      } else {
      oops::Log::info() << "Not doing adjoint test for inverse" << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVariableChangeInverse() {
  typedef LinearVariableChangeFixture<MODEL>       Test_;
  typedef oops::Increment<MODEL>                   Increment_;
  typedef oops::LinearVariableChangeBase<MODEL>    LinearVariableChange_;
  typedef oops::LinearVariableChangeFactory<MODEL> LinearVariableChangeFactory_;

  for (std::size_t jj = 0; jj < Test_::confs().size(); ++jj) {
    oops::Variables varin(Test_::confs()[jj], "input variables");
    oops::Variables varout(Test_::confs()[jj], "output variables");

    const double tol = Test_::confs()[jj].getDouble("tolerance inverse");

    const bool testinverse = Test_::confs()[jj].getBool("test inverse", false);
    if (testinverse)
      {
      oops::Log::info() << "Testing multiplyInverse" << std::endl;
      std::unique_ptr<LinearVariableChange_> changevar(LinearVariableChangeFactory_::create(
                                        Test_::xx(), Test_::xx(),
                                        Test_::resol(), Test_::confs()[jj]));

      Increment_    dxin(Test_::resol(), varin,  Test_::time());
      Increment_  KIdxin(Test_::resol(), varout, Test_::time());
      Increment_ KKIdxin(Test_::resol(), varin,  Test_::time());

      dxin.random();

      changevar->multiplyInverse(dxin, KIdxin);
      changevar->multiply(KIdxin, KKIdxin);

      const double zz1 = dxin.norm();
      const double zz2 = KKIdxin.norm();

      oops::Log::info() << "<x>, <KK^{-1}x>=" << zz1 << " " << zz2 << std::endl;
      oops::Log::info() << "<x>-<KK^{-1}x>=" << zz1-zz2 << std::endl;

      EXPECT((zz1-zz2) < tol);
    } else {
      oops::Log::info() << "multiplyInverse test not executed" << std::endl;
      EXPECT(1.0 < 2.0);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearVariableChange : public oops::Test {
 public:
  LinearVariableChange() {}
  virtual ~LinearVariableChange() {}
 private:
  std::string testid() const {return "test::LinearVariableChange<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/LinearVariableChange/testLinearVariableChangeZero")
      { testLinearVariableChangeZero<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearVariableChange/testLinearVariableChangeAdjoint")
      { testLinearVariableChangeAdjoint<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearVariableChange/testLinearVariableChangeInverse")
      { testLinearVariableChangeInverse<MODEL>(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LINEARVARIABLECHANGE_H_
