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

 private:
  static LinearVariableChangeFixture<MODEL>& getInstance() {
    static LinearVariableChangeFixture<MODEL> theLinearVariableChangeFixture;
    return theLinearVariableChangeFixture;
  }

  LinearVariableChangeFixture<MODEL>() {
    oops::instantiateVariableChangeFactory<MODEL>();

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
  typedef oops::LinearVariableChangeBase<MODEL>    LinearVariableChange_;
  typedef oops::LinearVariableChangeFactory<MODEL> LinearVariableChangeFactory_;

  for (std::size_t jj = 0; jj < Test_::confs().size(); ++jj) {
    oops::Variables varin(Test_::confs()[jj], "input variables");
    oops::Variables varout(Test_::confs()[jj], "output variables");

    std::unique_ptr<LinearVariableChange_> changevar(LinearVariableChangeFactory_::create(
                                      Test_::xx(), Test_::xx(),
                                      Test_::resol(), Test_::confs()[jj]));
    oops::Log::test() << "Testing linear variable change: " << *changevar << std::endl;
    Increment_  dxinTlIAd(Test_::resol(), varin,  Test_::time());
    Increment_  dxinAdInv(Test_::resol(), varout, Test_::time());
    Increment_ dxoutTlIAd(Test_::resol(), varout, Test_::time());
    Increment_ dxoutAdInv(Test_::resol(), varin,  Test_::time());

    // dxinTlIAd = 0, check if K.dxinTlIAd = 0
    dxinTlIAd.zero();
    dxoutTlIAd = changevar->multiply(dxinTlIAd);
    EXPECT(dxoutTlIAd.norm() == 0.0);

    // dxinAdInv = 0, check if K^T.dxinAdInv = 0
    dxinAdInv.zero();
    dxoutAdInv = changevar->multiplyAD(dxinAdInv);
    EXPECT(dxoutAdInv.norm() == 0.0);

    const bool testinverse = Test_::confs()[jj].getBool("test inverse", true);
    if (testinverse)
      {
        oops::Log::test() << "Doing zero test for inverse" << std::endl;
        dxinTlIAd.zero();
        dxoutTlIAd = changevar->multiplyInverseAD(dxinTlIAd);
        EXPECT(dxoutTlIAd.norm() == 0.0);

        dxinAdInv.zero();
        dxoutAdInv = changevar->multiplyInverse(dxinAdInv);
        EXPECT(dxoutAdInv.norm() == 0.0);
      } else {
      oops::Log::test() << "Not doing zero test for inverse" << std::endl;
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

    Increment_  dxinAdInv(Test_::resol(), varout, Test_::time());
    Increment_  dxinTlIAd(Test_::resol(), varin,  Test_::time());
    Increment_ dxoutAdInv(Test_::resol(), varin,  Test_::time());
    Increment_ dxoutTlIAd(Test_::resol(), varout, Test_::time());

    dxinAdInv.random();
    dxinTlIAd.random();

    Increment_  dxinAdInv0(dxinAdInv);
    Increment_  dxinTlIAd0(dxinTlIAd);

    dxoutTlIAd = changevar->multiply(dxinTlIAd);
    dxoutAdInv = changevar->multiplyAD(dxinAdInv);

    // zz1 = <dxoutTlIAd,dxinAdInv>
    double zz1 = dot_product(dxoutTlIAd, dxinAdInv0);
    // zz2 = <dxout,dxoutAdInv>
    double zz2 = dot_product(dxinTlIAd0, dxoutAdInv);

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
        dxoutAdInv.zero();
        dxoutTlIAd.zero();
        dxinAdInv.random();
        dxinTlIAd.random();
        dxinAdInv0 = dxinAdInv;
        dxinTlIAd0 = dxinTlIAd;
        dxoutTlIAd = changevar->multiplyInverseAD(dxinTlIAd);
        dxoutAdInv = changevar->multiplyInverse(dxinAdInv);
        zz1 = dot_product(dxoutTlIAd, dxinAdInv0);
        zz2 = dot_product(dxinTlIAd0, dxoutAdInv);
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
  typedef oops::LinearVariableChangeBase<MODEL>    LinearVariableChange_;
  typedef oops::LinearVariableChangeFactory<MODEL> LinearVariableChangeFactory_;

  for (std::size_t jj = 0; jj < Test_::confs().size(); ++jj) {
    oops::Variables varin(Test_::confs()[jj], "input variables");
    oops::Variables varout(Test_::confs()[jj], "output variables");

    const double tol = Test_::confs()[jj].getDouble("tolerance inverse");

    const bool testinverse = Test_::confs()[jj].getBool("test inverse", false);
    if (testinverse)
      {
      oops::Log::test() << "Testing multiplyInverse" << std::endl;
      std::unique_ptr<LinearVariableChange_> changevar(LinearVariableChangeFactory_::create(
                                        Test_::xx(), Test_::xx(),
                                        Test_::resol(), Test_::confs()[jj]));

      Increment_  dxinInv(Test_::resol(), varout, Test_::time());
      Increment_ dxoutInv(Test_::resol(), varin,  Test_::time());
      Increment_    dxout(Test_::resol(), varout, Test_::time());

      dxinInv.random();

      dxoutInv = changevar->multiplyInverse(dxinInv);
      dxout = changevar->multiply(dxoutInv);

      const double zz1 = dxinInv.norm();
      const double zz2 = dxout.norm();

      oops::Log::test() << "<x>, <KK^{-1}x>=" << zz1 << " " << zz2 << std::endl;
      oops::Log::test() << "<x>-<KK^{-1}x>=" << zz1-zz2 << std::endl;

      EXPECT((zz1-zz2) < tol);
    } else {
      oops::Log::test() << "multiplyInverse test not executed" << std::endl;
      EXPECT(1.0 < 2.0);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVariableChangeParametersWrapperValidName() {
  typedef LinearVariableChangeFixture<MODEL> Test_;
  for (const eckit::Configuration &config : Test_::confs()) {
    oops::LinearVariableChangeParametersWrapper<MODEL> parameters;
    EXPECT_NO_THROW(parameters.validateAndDeserialize(config));
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVariableChangeParametersWrapperInvalidName() {
  eckit::LocalConfiguration config;
  config.set("variable change", "###INVALID###");
  oops::LinearVariableChangeParametersWrapper<MODEL> parameters;
  if (oops::Parameters::isValidationSupported())
    EXPECT_THROWS_MSG(parameters.validate(config), "unrecognized enum value");
  EXPECT_THROWS_MSG(parameters.deserialize(config),
                    "does not exist in LinearVariableChangeFactory");
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearVariableChangeFactoryGetMakerNames() {
  typedef LinearVariableChangeFixture<MODEL> Test_;
  const std::vector<std::string> registeredNames =
      oops::LinearVariableChangeFactory<MODEL>::getMakerNames();
  for (const eckit::Configuration &config : Test_::confs()) {
    const std::string validName = config.getString("variable change");
    const bool found = std::find(registeredNames.begin(), registeredNames.end(), validName) !=
        registeredNames.end();
    EXPECT(found);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearVariableChange : public oops::Test {
 public:
  LinearVariableChange() {}
  virtual ~LinearVariableChange() {}

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
    ts.emplace_back(CASE("interface/LinearVariableChange/"
                         "testLinearVariableChangeParametersWrapperValidName")
      { testLinearVariableChangeParametersWrapperValidName<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearVariableChange/"
                         "testLinearVariableChangeParametersWrapperInvalidName")
      { testLinearVariableChangeParametersWrapperInvalidName<MODEL>(); });
    ts.emplace_back(CASE("interface/LinearVariableChange/"
                         "testLinearVariableChangeFactoryGetMakerNames")
      { testLinearVariableChangeFactoryGetMakerNames<MODEL>(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LINEARVARIABLECHANGE_H_
