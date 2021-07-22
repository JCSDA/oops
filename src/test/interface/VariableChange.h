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
#include "oops/generic/instantiateVariableChangeFactory.h"
#include "oops/interface/State.h"
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
    oops::instantiateVariableChangeFactory<MODEL>();

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
    VariableChange_ changevar(Test_::resol(), Test_::confs()[jj]);

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
    if (inverseFirst) {
      oops::Variables varin(Test_::confs()[jj], "input variables");
      State_ xin(Test_::resol(), varin, xx.validTime());
      changevar.changeVarInverse(xx, xin);
//      xx.zero();  Test for GEOS fails if uncommented
      changevar.changeVar(xin, xx);
    } else {
      oops::Variables varout(Test_::confs()[jj], "output variables");
      State_ xout(Test_::resol(), varout, xx.validTime());
      changevar.changeVar(xx, xout);
//      xx.zero();  Test for GEOS fails if uncommented
      changevar.changeVarInverse(xout, xx);
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

template <typename MODEL> void testVariableChangeParametersWrapperValidName() {
  typedef VariableChangeFixture<MODEL> Test_;
  for (const eckit::Configuration &config : Test_::confs()) {
    oops::VariableChangeParametersWrapper<MODEL> parameters;
    EXPECT_NO_THROW(parameters.validateAndDeserialize(config));
  }
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL> void testVariableChangeParametersWrapperInvalidName() {
  eckit::LocalConfiguration config;
  config.set("variable change", "###INVALID###");
  oops::VariableChangeParametersWrapper<MODEL> parameters;
  if (oops::Parameters::isValidationSupported())
    EXPECT_THROWS_MSG(parameters.validate(config), "unrecognized enum value");
  EXPECT_THROWS_MSG(parameters.deserialize(config),
                    "does not exist in VariableChangeFactory");
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL> void testVariableChangeFactoryGetMakerNames() {
  typedef VariableChangeFixture<MODEL> Test_;
  const std::vector<std::string> registeredNames =
      oops::VariableChangeFactory<MODEL>::getMakerNames();
  for (const eckit::Configuration &config : Test_::confs()) {
    const std::string validName = config.getString("variable change");
    const bool found = std::find(registeredNames.begin(), registeredNames.end(), validName) !=
        registeredNames.end();
    EXPECT(found);
  }
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL> class VariableChange : public oops::Test {
 public:
  VariableChange() {}
  virtual ~VariableChange() {VariableChangeFixture<MODEL>::reset();}
 private:
  std::string testid() const override {return "test::VariableChange<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/VariableChange/testVariableChangeInverse")
      { testVariableChangeInverse<MODEL>(); });
    ts.emplace_back(CASE("interface/VariableChange/"
                         "testVariableChangeParametersWrapperValidName")
      { testVariableChangeParametersWrapperValidName<MODEL>(); });
    ts.emplace_back(CASE("interface/VariableChange/"
                         "testVariableChangeParametersWrapperInvalidName")
      { testVariableChangeParametersWrapperInvalidName<MODEL>(); });
    ts.emplace_back(CASE("interface/VariableChange/"
                         "testVariableChangeFactoryGetMakerNames")
      { testVariableChangeFactoryGetMakerNames<MODEL>(); });
  }

  void clear() const override {}
};

// -------------------------------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_VARIABLECHANGE_H_
