/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_OBSOPERATOR_H_
#define TEST_INTERFACE_OBSOPERATOR_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

/// \brief tests constructor and print method
template <typename OBS> void testConstructor() {
  typedef oops::ObsOperator<OBS>             ObsOperator_;
  typedef ObsTestsFixture<OBS>               Test_;

  const std::vector<eckit::LocalConfiguration> obsConfs =
      TestEnvironment::config().getSubConfigurations("observations");

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::LocalConfiguration obsOpConf(obsConfs[jj], "obs operator");
    if (!obsConfs[jj].has("expect constructor to throw exception with message")) {
      auto hop = std::make_unique<ObsOperator_>(Test_::obspace()[jj], obsOpConf);
      EXPECT(hop.get());
      oops::Log::info() << "Testing ObsOperator: " << *hop << std::endl;
      hop.reset();
      EXPECT(!hop.get());
    } else {
      // The constructor is expected to throw an exception containing the specified string.
      const std::string expectedMessage =
          obsConfs[jj].getString("expect constructor to throw exception with message");
      EXPECT_THROWS_MSG(ObsOperator_(Test_::obspace()[jj], obsOpConf),
                        expectedMessage.c_str());
    }
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testSimulateObs() {
  typedef oops::GeoVaLs<OBS>           GeoVaLs_;
  typedef oops::ObsDiagnostics<OBS>    ObsDiags_;
  typedef oops::ObsAuxControl<OBS>     ObsAuxCtrl_;
  typedef oops::ObsOperator<OBS>       ObsOperator_;
  typedef oops::ObsVector<OBS>         ObsVector_;
  typedef ObsTestsFixture<OBS>         Test_;

  const std::vector<eckit::LocalConfiguration> obsConfs =
      TestEnvironment::config().getSubConfigurations("observations");

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::Configuration & obsConf = obsConfs[jj];
    if (obsConf.has("expect constructor to throw exception with message"))
      continue;

    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)

    oops::ObsDataVector<OBS, int> qc_flags(
      Test_::obspace()[jj],
      Test_::obspace()[jj].obsvariables(),
      std::string());
    const eckit::LocalConfiguration obsOpConf(obsConf, "obs operator");
    ObsOperator_ hop(Test_::obspace()[jj], obsOpConf);

    // initialize bias correction
    const eckit::LocalConfiguration bconf = obsConf.getSubConfiguration("obs bias");
    ObsAuxCtrl_ ybias(Test_::obspace()[jj], bconf);

    // initialize geovals
    oops::Variables hopvars = hop.requiredVars();
    oops::Variables reducedHopvars = ybias.requiredVars();
    hopvars += reducedHopvars;  // the reduced format is derived from the sampled format
    // read geovals from the file (in the sampled format)
    GeoVaLs_ gval(eckit::LocalConfiguration(obsConf, "geovals"), Test_::obspace()[jj], hopvars);
    // convert geovals to the reduced format
    hop.computeReducedVars(reducedHopvars, gval);

    // create obsvector to hold H(x)
    ObsVector_ hofx(Test_::obspace()[jj]);

    // create obsvector to hold bias
    ObsVector_ bias(Test_::obspace()[jj]);
    bias.zero();

    // create diagnostics to hold HofX diags
    oops::ObsVariables diagvars;
    diagvars += ybias.requiredHdiagnostics();
    ObsDiags_ diags(Test_::obspace()[jj], hop.locations(), diagvars);

    // call H(x), save result in the output file as @hofx
    if (obsConf.has("expect simulateObs to throw exception with message")) {
      // The simulateObs method is expected to throw an exception
      // containing the specified string.
      const std::string expectedMessage =
          obsConf.getString("expect simulateObs to throw exception with message");
      EXPECT_THROWS_MSG(hop.simulateObs(gval, hofx, ybias, qc_flags, bias, diags),
                        expectedMessage.c_str());
      continue;
    } else {
      hop.simulateObs(gval, hofx, ybias, qc_flags, bias, diags);
    }
    hofx.save("hofx");
    bias.save("ObsBias");

    const double tol = obsConf.getDouble("tolerance");
    if (obsConf.has("vector ref")) {
      // if reference h(x) is saved in file as a vector, read from file
      // and compare the norm of difference to zero
      ObsVector_ obsref(Test_::obspace()[jj], obsConf.getString("vector ref"));
      obsref -= hofx;
      const double zz = obsref.rms();
      oops::Log::info() << "Vector difference between reference and computed: " << obsref;
      EXPECT(zz < 100*tol);  //  change tol from percent to actual value.
                             //  tol used in is_close is relative
    } else if (obsConf.has("norm ref")) {
      // if reference h(x) is saved in file as a vector, read from file
      // and compare the difference, normalised by the reference values to zero
      ObsVector_ obsref(Test_::obspace()[jj], obsConf.getString("norm ref"));
      obsref -= hofx;
      obsref /= hofx;
      const double zz = obsref.rms();
      oops::Log::info() << "Normalised vector difference between reference and computed: "
                        << obsref;
      EXPECT(zz < 100*tol);  //  change tol from percent to actual value.
                             //  tol used in is_close is relative
    } else {
      // else compare h(x) norm to the norm from the config
      const double zz = hofx.rms();
      const double xx = obsConf.getDouble("rms ref");
      EXPECT(oops::is_close(xx, zz, tol));
    }
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsOperator : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;
 public:
  using oops::Test::Test;
  virtual ~ObsOperator() {}
 private:
  std::string testid() const override {return "test::ObsOperator<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsOperator/testConstructor")
      { testConstructor<OBS>(); });
    ts.emplace_back(CASE("interface/ObsOperator/testSimulateObs")
      { testSimulateObs<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSOPERATOR_H_
