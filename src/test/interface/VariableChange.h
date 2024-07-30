/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_VARIABLECHANGE_H_
#define TEST_INTERFACE_VARIABLECHANGE_H_

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
#include "oops/base/State.h"
#include "oops/base/VariableChangeParametersBase.h"
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

namespace test {

// -------------------------------------------------------------------------------------------------

template <typename MODEL> class VariableChangeFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>  Geometry_;
  typedef oops::State<MODEL>     State_;

 public:
  static std::vector<eckit::LocalConfiguration> & confs() {return getInstance().confs_;}
  static const Geometry_       & resol()  {return *getInstance().resol_;}
  static void reset() {
    getInstance().resol_.reset();
  }

 private:
  static VariableChangeFixture<MODEL>& getInstance() {
    static VariableChangeFixture<MODEL> theVariableChangeFixture;
    return theVariableChangeFixture;
  }

  VariableChangeFixture<MODEL>() {
    // Geometry for the test
    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::world()));

    // Configuration (list of all variable changes)
    TestEnvironment::config().get("variable change tests", confs_);
  }

  ~VariableChangeFixture<MODEL>() {}

  std::vector<eckit::LocalConfiguration>  confs_;
  std::unique_ptr<const Geometry_>        resol_;
};

// -------------------------------------------------------------------------------------------------

template <typename MODEL> void testVariableChangeInverse() {
  typedef VariableChangeFixture<MODEL>   Test_;
  typedef oops::State<MODEL>             State_;
  typedef oops::VariableChange<MODEL>    VariableChange_;

  // Loop over all variable changes
  for (std::size_t jj = 0; jj < Test_::confs().size(); ++jj) {
    // Construct variable change
    const eckit::LocalConfiguration changeVarConfig(Test_::confs()[jj], "variable change");
    VariableChange_ changevar(changeVarConfig, Test_::resol());

    oops::Log::test() << "Testing VariableChange: " << changevar << std::endl;
    // User specified tolerance for pass/fail
    const double tol = Test_::confs()[jj].getDouble("tolerance inverse");

    // Create states with input and output variables
    const eckit::LocalConfiguration initialConfig(Test_::confs()[jj], "state");
    State_ xx(Test_::resol(), initialConfig);

    const double xxnorm_ref = xx.norm();

    // Order, inverse first or not (default)
    // Note: switch input and output variables in configuration if true
    const bool inverseFirst = Test_::confs()[jj].getBool("inverse first", false);

    // Convert from input to output variables and back (or vice versa)
    oops::Variables varin(changeVarConfig, "input variables");
    oops::Variables varout(changeVarConfig, "output variables");
    if (inverseFirst) {
      changevar.changeVarInverse(xx, varin);
      changevar.changeVar(xx, varout);
    } else {
      changevar.changeVar(xx, varout);
      oops::Log::debug() << "Test output of changeVar: " << xx << std::endl;
      changevar.changeVarInverse(xx, varin);
      oops::Log::debug() << "Test output of changeVarInverse: " << xx << std::endl;
    }

    // Compute norms of the result and reference
    const double xxnorm_tst =   xx.norm();

    // Print the input and final state
    oops::Log::test() << "<xin>, <K^{-1}[K(xin)]>, (<xin>-<K^{-1}[K(xin)]<xin>)/>=" << xxnorm_ref <<
                      " " << xxnorm_tst << " " << (xxnorm_ref - xxnorm_tst)/xxnorm_ref <<std::endl;

    // Is result similar to the reference
    EXPECT(oops::is_close(xxnorm_tst, xxnorm_ref, tol));
  }
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL> void testVariableChange() {
  typedef VariableChangeFixture<MODEL>   Test_;
  typedef oops::State<MODEL>             State_;
  typedef oops::VariableChange<MODEL>    VariableChange_;

  // Loop over all variable changes
  for (std::size_t jj = 0; jj < Test_::confs().size(); ++jj) {
    // Construct variable change
    const eckit::LocalConfiguration changeVarConfig(Test_::confs()[jj], "variable change");
    VariableChange_ changevar(changeVarConfig, Test_::resol());

    oops::Log::test() << "Testing VariableChange: " << changevar << std::endl;

    // Create state with input variables
    const eckit::LocalConfiguration initialConfig(Test_::confs()[jj], "state");
    State_ xx(Test_::resol(), initialConfig);
    oops::Log::test() << "State before variable change: " << xx << std::endl;

    // Change variables and print output state
    oops::Variables varout(changeVarConfig, "output variables");
    oops::Log::test() << "Change of variables to " << varout.variables() << std::endl;
    changevar.changeVar(xx, varout);
    oops::Log::test() << "State after variable change: " << xx << std::endl;
  }
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL> void testVariableChangeParametersValidName() {
  typedef VariableChangeFixture<MODEL> Test_;
  typedef oops::VariableChange<MODEL>  VariableChange_;
  typedef typename VariableChange_::Parameters_ Parameters_;
  for (const eckit::Configuration &config : Test_::confs()) {
    Parameters_ parameters;
    const eckit::LocalConfiguration changeVarConfig(config, "variable change");
    EXPECT_NO_THROW(parameters.validateAndDeserialize(changeVarConfig));
  }
}

// -------------------------------------------------------------------------------------------------
template <typename MODEL> class VariableChange : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~VariableChange() {VariableChangeFixture<MODEL>::reset();}
 private:
  std::string testid() const override {return "test::VariableChange<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    // The testVariableChangeInverse test is broken for some models due to the recent
    // interface change to the changeVar method. (The old interface could work
    // work with non-invertable variables, but the new one cannot.) The test is
    // removed until a better test can be written.
    //
    // ts.emplace_back(CASE("interface/VariableChange/testVariableChangeInverse")
    //   { testVariableChangeInverse<MODEL>(); });
    ts.emplace_back(CASE("interface/VariableChange/testVariableChange")
      { testVariableChange<MODEL>(); });
    ts.emplace_back(CASE("interface/VariableChange/testVariableChangeParametersValidName")
      { testVariableChangeParametersValidName<MODEL>(); });
  }

  void clear() const override {}
};

// -------------------------------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_VARIABLECHANGE_H_
