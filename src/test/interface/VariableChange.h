/*
 * (C) Copyright 2018  UCAR
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
#include "oops/base/VariableChangeBase.h"
#include "oops/generic/instantiateVariableChangeFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

namespace test {

// -------------------------------------------------------------------------------------------------

template <typename MODEL> class VariableChangeFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>  Geometry_;
  typedef oops::State<MODEL>     State_;

 public:
  static std::vector<eckit::LocalConfiguration> & confs() {return getInstance().confs_;}
  static const State_          & xx()     {return *getInstance().xx_;}
  static const Geometry_       & resol()  {return *getInstance().resol_;}

 private:
  static VariableChangeFixture<MODEL>& getInstance() {
    static VariableChangeFixture<MODEL> theVariableChangeFixture;
    return theVariableChangeFixture;
  }

  VariableChangeFixture<MODEL>() {
    oops::instantiateVariableChangeFactory<MODEL>();

    // Geometry for the test
    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::comm()));

    // Configuration (list of all variable changes)
    TestEnvironment::config().get("VariableChangeTests", confs_);
  }

  ~VariableChangeFixture<MODEL>() {}

  std::vector<eckit::LocalConfiguration>  confs_;
  std::unique_ptr<const Geometry_>        resol_;
};

// -------------------------------------------------------------------------------------------------

template <typename MODEL> void testVariableChangeInverse() {
  typedef VariableChangeFixture<MODEL>   Test_;
  typedef oops::State<MODEL>                   State_;
  typedef oops::VariableChangeBase<MODEL>      VariableChange_;
  typedef oops::VariableChangeFactory<MODEL>   VariableChangeFactory_;

  // Loop over all variable changes
  for (std::size_t jj = 0; jj < Test_::confs().size(); ++jj) {
    eckit::LocalConfiguration varinconf(Test_::confs()[jj], "inputVariables");
    eckit::LocalConfiguration varoutconf(Test_::confs()[jj], "outputVariables");
    oops::Variables varin(varinconf);
    oops::Variables varout(varoutconf);

    // Construct variable change
    std::unique_ptr<VariableChange_> \
      changevar(VariableChangeFactory_::create(Test_::confs()[jj], Test_::resol()));

    // User specified tolerance for pass/fail
    const double tol = Test_::confs()[jj].getDouble("toleranceInverse");

    // Create states with input and output variables
    const eckit::LocalConfiguration initialConfig(Test_::confs()[jj], "state");
    State_  xin(Test_::resol(), varin, initialConfig);
    State_ xout(Test_::resol(), varout, xin.validTime());

    // Save copy of the initial state
    State_ xref(xin);

    // Order, inverse first or not (default)
    // Note: switch input and output variables in configuration if true
    const bool inverseFirst = Test_::confs()[jj].getBool("inverseFirst", false);

    // Convert from input to output variables and back (or vice versa)
    if (inverseFirst) {
      changevar->changeVarInverse(xin, xout);
      changevar->changeVar(xout, xin);
    } else {
      changevar->changeVar(xin, xout);
      changevar->changeVarInverse(xout, xin);
    }

    // Compute norms of the result and reference
    const double xxnorm_ref = xref.norm();
    const double xxnorm_tst =  xin.norm();

    // Print the input and final state
    oops::Log::info() << "<xin>, <K^{-1}[K(xin)]>, (<xin>-<K^{-1}[K(xin)]<xin>)/>=" << xxnorm_ref <<
                      " " << xxnorm_tst << " " << (xxnorm_ref - xxnorm_tst)/xxnorm_ref <<std::endl;

    // Is result similar to the reference
    EXPECT(oops::is_close(xxnorm_tst, xxnorm_ref, tol));
  }
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL> class VariableChange : public oops::Test {
 public:
  VariableChange() {}
  virtual ~VariableChange() {}
 private:
  std::string testid() const {return "test::VariableChange<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/VariableChange/testVariableChangeInverse")
      { testVariableChangeInverse<MODEL>(); });
  }
};

// -------------------------------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_VARIABLECHANGE_H_
