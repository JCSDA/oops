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
    oops::Log::info() << "Testing linear variable change" << std::endl;
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

    oops::Log::info() << "<dxout,KTdxin>-<Kdxout,dxin>/<dxout,KTdxin>="
                      << (zz1-zz2)/zz1 << std::endl;
    oops::Log::info() << "<dxout,KTdxin>-<Kdxout,dxin>/<Kdxout,dxin>="
                      << (zz1-zz2)/zz2 << std::endl;
    const double tol = 1e-10;
    EXPECT(oops::is_close(zz1, zz2, tol));
  }
}

// -----------------------------------------------------------------------------

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
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LINEARVARIABLECHANGE_H_
