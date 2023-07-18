/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_ERRORCOVARIANCE_H_
#define TEST_INTERFACE_ERRORCOVARIANCE_H_

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
#include "oops/base/Increment4D.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

namespace test {

// =============================================================================

template <typename MODEL> class ErrorCovarianceFixture : private boost::noncopyable {
  typedef oops::ModelSpaceCovarianceBase<MODEL> Covariance_;
  typedef oops::Geometry<MODEL>       Geometry_;

 public:
  static const eckit::Configuration & test()       {return *getInstance().test_;}
  static const Geometry_            & resol()      {return *getInstance().resol_;}
  static const oops::Variables      & ctlvars()    {return *getInstance().ctlvars_;}
  static const std::vector<util::DateTime> & time() {return getInstance().time_;}
  static const Covariance_          & covariance() {return *getInstance().B_;}
  static void reset() {
    getInstance().B_.reset();
    getInstance().time_.clear();
    getInstance().ctlvars_.reset();
    getInstance().resol_.reset();
    getInstance().test_.reset();
  }

 private:
  static ErrorCovarianceFixture<MODEL>& getInstance() {
    static ErrorCovarianceFixture<MODEL> theErrorCovarianceFixture;
    return theErrorCovarianceFixture;
  }

  ErrorCovarianceFixture<MODEL>() {
    oops::instantiateCovarFactory<MODEL>();

    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "covariance test"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::world()));

    ctlvars_.reset(new oops::Variables(TestEnvironment::config(), "analysis variables"));

    const eckit::LocalConfiguration fgconf(TestEnvironment::config(), "background");
    oops::State4D<MODEL> xx(*resol_, fgconf);

    time_ = xx.times();

//  Setup the B matrix
    const eckit::LocalConfiguration covar(TestEnvironment::config(), "background error");
    B_.reset(oops::CovarianceFactory<MODEL>::create(*resol_, *ctlvars_, covar, xx, xx));
  }

  ~ErrorCovarianceFixture<MODEL>() {}

  std::unique_ptr<const eckit::LocalConfiguration>   test_;
  std::unique_ptr<const Geometry_>       resol_;
  std::unique_ptr<const oops::Variables> ctlvars_;
  std::vector<util::DateTime>            time_;
  std::unique_ptr<Covariance_>           B_;
};

// =============================================================================

template <typename MODEL> void testErrorCovarianceZero() {
  typedef ErrorCovarianceFixture<MODEL>   Test_;
  typedef oops::Increment4D<MODEL>    Increment_;

  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  Increment_ dx2(Test_::resol(), Test_::ctlvars(), Test_::time());

  Test_::covariance().randomize(dx2);
  oops::Log::test() << "dx2.norm()=" << dx2[0].norm() << std::endl;

  EXPECT(dx1[0].norm() == 0.0);
  EXPECT(dx2[0].norm() > 0.0);

  Test_::covariance().multiply(dx1, dx2);
  EXPECT(dx2[0].norm() == 0.0);

  const bool testinverse = Test_::test().getBool("testinverse", true);
  if (testinverse)
    {
      oops::Log::test() << "Doing zero test for inverse" << std::endl;
      dx1.zero();
      Test_::covariance().randomize(dx2);
      EXPECT(dx1[0].norm() == 0.0);
      EXPECT(dx2[0].norm() > 0.0);
      Test_::covariance().inverseMultiply(dx1, dx2);
      EXPECT(dx2[0].norm() == 0.0);
    } else {
      oops::Log::test() << "Not doing zero test for inverse" << std::endl;
    }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testErrorCovarianceInverse() {
  typedef ErrorCovarianceFixture<MODEL>   Test_;
  typedef oops::Increment4D<MODEL>    Increment_;

  const bool testinverse = Test_::test().getBool("testinverse", true);
  if (testinverse)
    {
      Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
      Increment_ dx2(Test_::resol(), Test_::ctlvars(), Test_::time());
      Increment_ dx3(Test_::resol(), Test_::ctlvars(), Test_::time());
      Test_::covariance().randomize(dx1);
      EXPECT(dx1[0].norm() > 0.0);

      Test_::covariance().multiply(dx1, dx2);
      Test_::covariance().inverseMultiply(dx2, dx3);

      EXPECT(dx2[0].norm() > 0.0);
      EXPECT(dx3[0].norm() > 0.0);
      dx3 -= dx1;
      const double tol = Test_::test().getDouble("tolerance");
      EXPECT(dx3[0].norm()/dx1[0].norm() < tol);
    } else {
      return;
    }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testErrorCovarianceSym() {
  typedef ErrorCovarianceFixture<MODEL>   Test_;
  typedef oops::Increment4D<MODEL>    Increment_;

  Increment_ dx(Test_::resol(), Test_::ctlvars(), Test_::time());
  Increment_ Bdx(Test_::resol(), Test_::ctlvars(), Test_::time());
  Increment_ dy(Test_::resol(), Test_::ctlvars(), Test_::time());
  Increment_ Bdy(Test_::resol(), Test_::ctlvars(), Test_::time());

  dx.random();
  dy.random();

  Test_::covariance().multiply(dx, Bdx);
  Test_::covariance().multiply(dy, Bdy);
  const double zz1 = dot_product(dx, Bdy);
  const double zz2 = dot_product(Bdx, dy);
  oops::Log::test() << "<dx,Bdy>-<Bdx,dy>/<dx,Bdy>="
                    <<  (zz1-zz2)/zz1 << std::endl;
  oops::Log::test() << "<dx,Bdy>-<Bdx,dy>/<Bdx,dy>="
                    <<  (zz1-zz2)/zz2 << std::endl;
  const double tol = Test_::test().getDouble("tolerance");
  EXPECT(oops::is_close(zz1, zz2, tol));
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCovarianceParametersWrapperValidName() {
  eckit::LocalConfiguration config(TestEnvironment::config(), "background error");
  oops::ModelSpaceCovarianceParametersWrapper<MODEL> parameters;
  EXPECT_NO_THROW(parameters.validateAndDeserialize(config));
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCovarianceParametersWrapperInvalidName() {
  eckit::LocalConfiguration config;
  config.set("covariance model", "###INVALID###");
  oops::ModelSpaceCovarianceParametersWrapper<MODEL> parameters;
  if (oops::Parameters::isValidationSupported())
    EXPECT_THROWS_MSG(parameters.validate(config), "unrecognized enum value");
  EXPECT_THROWS_MSG(parameters.deserialize(config),
                    "does not exist in CovarianceFactory");
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCovarianceFactoryGetMakerNames() {
  eckit::LocalConfiguration config(TestEnvironment::config(), "background error");
  const std::string validName = config.getString("covariance model");
  const std::vector<std::string> registeredNames =
      oops::CovarianceFactory<MODEL>::getMakerNames();
  const bool found = std::find(registeredNames.begin(), registeredNames.end(), validName) !=
      registeredNames.end();
  EXPECT(found);
}

// =============================================================================

template <typename MODEL>
class ErrorCovariance : public oops::Test  {
 public:
  ErrorCovariance() {}
  virtual ~ErrorCovariance() {ErrorCovarianceFixture<MODEL>::reset();}

 private:
  std::string testid() const override {return "test::ErrorCovariance<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ErrorCovariance/testErrorCovarianceZero")
      { testErrorCovarianceZero<MODEL>(); });
    ts.emplace_back(CASE("interface/ErrorCovariance/testErrorCovarianceInverse")
      { testErrorCovarianceInverse<MODEL>(); });
    ts.emplace_back(CASE("interface/ErrorCovariance/testErrorCovarianceSym")
      { testErrorCovarianceSym<MODEL>(); });
    ts.emplace_back(CASE("interface/ErrorCovariance/testCovarianceParametersWrapperValidName")
      { testCovarianceParametersWrapperValidName<MODEL>(); });
    ts.emplace_back(CASE("interface/ErrorCovariance/testCovarianceParametersWrapperInvalidName")
      { testCovarianceParametersWrapperInvalidName<MODEL>(); });
    ts.emplace_back(CASE("interface/ErrorCovariance/testCovarianceFactoryGetMakerNames")
      { testCovarianceFactoryGetMakerNames<MODEL>(); });
  }

  void clear() const override {}
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_ERRORCOVARIANCE_H_
