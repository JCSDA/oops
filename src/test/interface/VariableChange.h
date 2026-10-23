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

template <typename MODEL> class VariableChange : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~VariableChange() {VariableChangeFixture<MODEL>::reset();}
 private:
  std::string testid() const override {return "test::VariableChange<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/VariableChange/testVariableChange")
      { testVariableChange<MODEL>(); });
  }

  void clear() const override {}
};

// -------------------------------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_VARIABLECHANGE_H_
