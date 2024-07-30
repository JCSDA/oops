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
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Expect.h"
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
  static void reset() {
    getInstance().time_.reset();
    getInstance().xx_.reset();
    getInstance().resol_.reset();
  }

 private:
  static LinearVariableChangeFixture<MODEL>& getInstance() {
    static LinearVariableChangeFixture<MODEL> theLinearVariableChangeFixture;
    return theLinearVariableChangeFixture;
  }

  LinearVariableChangeFixture<MODEL>() {
    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::world()));

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
  typedef oops::LinearVariableChange<MODEL>    LinearVariableChange_;

  for (std::size_t jj = 0; jj < Test_::confs().size(); ++jj) {
    const eckit::LocalConfiguration lvcConfig(Test_::confs()[jj], "linear variable change");
    oops::Variables varin(lvcConfig, "input variables");
    oops::Variables varout(lvcConfig, "output variables");

    std::unique_ptr<LinearVariableChange_>
        changevar(new LinearVariableChange_(Test_::resol(), lvcConfig));
    oops::Log::test() << "Testing linear variable change" << std::endl;
    Increment_  dxChangeVarTL(Test_::resol(), varin,  Test_::time());
    Increment_  dxChangeVarAD(Test_::resol(), varout,  Test_::time());

    changevar->changeVarTraj(Test_::xx(), varout);

    // dxMultiply = 0, check if K.dxMultiply = 0
    dxChangeVarTL.zero();
    changevar->changeVarTL(dxChangeVarTL, varout);
    EXPECT(dxChangeVarTL.norm() == 0.0);

    // dxinAdInv = 0, check if K^T.dxinAdInv = 0
    dxChangeVarAD.zero();
    changevar->changeVarAD(dxChangeVarAD, varin);
    EXPECT(dxChangeVarAD.norm() == 0.0);

    const bool testinverse = Test_::confs()[jj].getBool("test inverse", true);
    if (testinverse)
      {
        oops::Log::test() << "Doing zero test for inverse" << std::endl;
        Increment_  dxChangeVarInverseTL(Test_::resol(), varout,  Test_::time());
        Increment_  dxChangeVarInverseAD(Test_::resol(), varin,  Test_::time());
        dxChangeVarInverseAD.zero();
        changevar->changeVarInverseAD(dxChangeVarInverseAD, varout);
        EXPECT(dxChangeVarInverseAD.norm() == 0.0);

        dxChangeVarInverseTL.zero();
        changevar->changeVarInverseTL(dxChangeVarInverseTL, varin);
        EXPECT(dxChangeVarInverseTL.norm() == 0.0);
      } else {
      oops::Log::test() << "Not doing zero test for inverse" << std::endl;
    }
  }
}
// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVariableChangeAdjoint() {
  typedef LinearVariableChangeFixture<MODEL>       Test_;
  typedef oops::Increment<MODEL>                   Increment_;
  typedef oops::LinearVariableChange<MODEL>    LinearVariableChange_;

  for (std::size_t jj = 0; jj < Test_::confs().size(); ++jj) {
    const eckit::LocalConfiguration lvcConfig(Test_::confs()[jj], "linear variable change");
    oops::Variables varin(lvcConfig, "input variables");
    oops::Variables varout(lvcConfig, "output variables");

    std::unique_ptr<LinearVariableChange_>
        changevar(new LinearVariableChange_(Test_::resol(), lvcConfig));
    changevar->changeVarTraj(Test_::xx(), varout);

    Increment_  dxChangeVarTL(Test_::resol(), varin,  Test_::time());
    Increment_  dxChangeVarAD(Test_::resol(), varout,  Test_::time());

    dxChangeVarTL.random();
    dxChangeVarAD.random();

    Increment_  dxChangeVarTLIn(dxChangeVarTL);
    Increment_  dxChangeVarADIn(dxChangeVarAD);

    changevar->changeVarTL(dxChangeVarTL, varout);
    changevar->changeVarAD(dxChangeVarAD, varin);

    double zz1 = dot_product(dxChangeVarTL, dxChangeVarADIn);
    double zz2 = dot_product(dxChangeVarTLIn, dxChangeVarAD);

    oops::Log::test() << "<dxout,KTdxin>-<Kdxout,dxin>/<dxout,KTdxin>="
                      << (zz1-zz2)/zz1 << std::endl;
    oops::Log::test() << "<dxout,KTdxin>-<Kdxout,dxin>/<Kdxout,dxin>="
                      << (zz1-zz2)/zz2 << std::endl;
    const double tol = 1e-10;
    EXPECT(oops::is_close(zz1, zz2, tol));
    const bool testinverse = Test_::confs()[jj].getBool("test inverse", true);
    if (testinverse)
    {
      oops::Log::test() << "Doing adjoint test for inverse" << std::endl;
      Increment_  dxChangeVarInverseTL(Test_::resol(), varout,  Test_::time());
      Increment_  dxChangeVarInverseAD(Test_::resol(), varin,  Test_::time());
      dxChangeVarInverseTL.random();
      dxChangeVarInverseAD.random();

      Increment_  dxChangeVarInverseTLIn(dxChangeVarInverseTL);
      Increment_  dxChangeVarInverseADIn(dxChangeVarInverseAD);

      changevar->changeVarInverseAD(dxChangeVarInverseAD, varout);
      changevar->changeVarInverseTL(dxChangeVarInverseTL, varin);

      zz1 = dot_product(dxChangeVarInverseTL, dxChangeVarInverseADIn);
      zz2 = dot_product(dxChangeVarInverseTLIn, dxChangeVarInverseAD);
      oops::Log::test() << "<dxout,KinvTdxin>-<Kinvdxout,dxin>/<dxout,KinvTdxin>="
                      << (zz1-zz2)/zz1 << std::endl;
      oops::Log::test() << "<dxout,KinvTdxin>-<Kinvdxout,dxin>/<Kinvdxout,dxin>="
                      << (zz1-zz2)/zz2 << std::endl;
      EXPECT(oops::is_close(zz1, zz2, tol));
    } else {
      oops::Log::test() << "Not doing adjoint test for inverse" << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVariableChangeInverse() {
  typedef LinearVariableChangeFixture<MODEL>       Test_;
  typedef oops::Increment<MODEL>                   Increment_;
  typedef oops::LinearVariableChange<MODEL>    LinearVariableChange_;

  for (std::size_t jj = 0; jj < Test_::confs().size(); ++jj) {
    const eckit::LocalConfiguration lvcConfig(Test_::confs()[jj], "linear variable change");
    oops::Variables varin(lvcConfig, "input variables");
    oops::Variables varout(lvcConfig, "output variables");

    const double tol = Test_::confs()[jj].getDouble("tolerance inverse");

    const bool testinverse = Test_::confs()[jj].getBool("test inverse", false);
    if (testinverse)
    {
      oops::Log::test() << "Testing multiplyInverse" << std::endl;
      std::unique_ptr<LinearVariableChange_>
        changevar(new LinearVariableChange_(Test_::resol(), lvcConfig));
      changevar->changeVarTraj(Test_::xx(), varout);

      Increment_  dx(Test_::resol(), varout, Test_::time());
      dx.random();
      Increment_  dx0(dx);

      changevar->changeVarInverseTL(dx, varin);
      changevar->changeVarTL(dx, varout);

      const double zz1 = dx.norm();
      const double zz2 = dx0.norm();

      oops::Log::test() << "<x>, <KK^{-1}x>=" << zz1 << " " << zz2 << std::endl;
      oops::Log::test() << "<x>-<KK^{-1}x>=" << zz1-zz2 << std::endl;

      EXPECT((zz1-zz2) < tol);
    } else {
      oops::Log::test() << "multiplyInverse test not executed" << std::endl;
      EXPECT(true);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVarChangeParametersValidName() {
  typedef LinearVariableChangeFixture<MODEL> Test_;
  typedef oops::LinearVariableChange<MODEL>  LinearVariableChange_;
  for (const eckit::Configuration &config : Test_::confs()) {
    typename LinearVariableChange_::Parameters_ parameters;
    const eckit::LocalConfiguration lvcConfig(config, "linear variable change");
    EXPECT_NO_THROW(parameters.validateAndDeserialize(lvcConfig));
  }
}

// -------------------------------------------------------------------------------------------------


template <typename MODEL>
class LinearVariableChange : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~LinearVariableChange() {LinearVariableChangeFixture<MODEL>::reset();}

 private:
  std::string testid() const override {return "test::LinearVariableChange<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/LinearVariableChange/testLinearVariableChangeZero")
      { testLinearVariableChangeZero<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearVariableChange/testLinearVariableChangeAdjoint")
      { testLinearVariableChangeAdjoint<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearVariableChange/testLinearVariableChangeInverse")
      { testLinearVariableChangeInverse<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearVariableChange/testLinearVarChangeParametersValidName")
      { testLinearVarChangeParametersValidName<MODEL>(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LINEARVARIABLECHANGE_H_
