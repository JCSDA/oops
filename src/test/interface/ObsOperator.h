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
#include "oops/base/ObsTypeParameters.h"
#include "oops/base/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

/// \brief Options used to configure a test simulating observations from a single obs space
/// using a particular ObsOperator.
template <typename OBS>
class ObsTypeParameters : public oops::ObsTypeParametersBase<OBS> {
  typedef oops::ObsTypeParametersBase<OBS> Base;
  OOPS_CONCRETE_PARAMETERS(ObsTypeParameters, Base)

  typedef typename oops::GeoVaLs<OBS>::Parameters_ GeoVaLsParameters_;

 public:
  /// Overridden to perform some extra checks that can't currently be done by the JSON schema
  /// validator.
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

 public:
  /// Options used to load GeoVaLs from a file.
  oops::RequiredParameter<GeoVaLsParameters_> geovals{"geovals", this};

  // One of these parameters must be set.
  oops::OptionalParameter<std::string> expectConstructorToThrow{
    "expect constructor to throw exception with message", this};
  oops::OptionalParameter<std::string> expectSimulateObsToThrow{
    "expect simulateObs to throw exception with message", this};
  oops::OptionalParameter<double> tolerance{"tolerance", this};

  // One of these parameters must be set if `tolerance` is.
  oops::OptionalParameter<std::string> vectorRef{"vector ref", this};
  oops::OptionalParameter<std::string> normRef{"norm ref", this};
  oops::OptionalParameter<double> rmsRef{"rms ref", this};

 private:
  // Parameters ignored by this test but used by the LinearObsOperator test. Both tests tend to
  // use the same YAML files.
  oops::Parameter<eckit::LocalConfiguration> linearObsOperatorTest{
    "linear obs operator test", eckit::LocalConfiguration(), this};
  oops::Parameter<eckit::LocalConfiguration> expectSetTrajectoryToThrow{
    "expect setTrajectory to throw exception with message", eckit::LocalConfiguration(), this};
  oops::Parameter<eckit::LocalConfiguration> expectSimulateObsTLToThrow{
    "expect simulateObsTL to throw exception with message", eckit::LocalConfiguration(), this};
  oops::Parameter<eckit::LocalConfiguration> expectSimulateObsADToThrow{
    "expect simulateObsAD to throw exception with message", eckit::LocalConfiguration(), this};
};

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsTypeParameters<OBS>::deserialize(util::CompositePath &path,
                                         const eckit::Configuration &config) {
  Base::deserialize(path, config);

  if (tolerance.value() == boost::none &&
      expectConstructorToThrow.value() == boost::none &&
      expectSimulateObsToThrow.value() == boost::none) {
    throw eckit::UserError(path.path() +
                           ": At least one of the 'tolerance', "
                           "'expect constructor to throw exception with message' and "
                           "'expect simulateObs to throw exception with message' options "
                           "must be set");
  }

  if (tolerance.value() != boost::none && vectorRef.value() == boost::none &&
      normRef.value() == boost::none && rmsRef.value() == boost::none) {
    throw eckit::UserError(path.path() +
                           ": If the 'tolerance' option is set, one of the 'vector ref', "
                           "'norm ref' and 'rms ref' options must be set as well");
  }
}

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the ObsOperator test.
template <typename OBS>
class TestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TestParameters, Parameters)

 public:
  /// Only observations taken at times lying in the (`window begin`, `window end`] interval
  /// will be included in observation spaces.
  oops::RequiredParameter<util::DateTime> windowBegin{"window begin", this};
  oops::RequiredParameter<util::DateTime> windowEnd{"window end", this};
  /// Each element of this list configures an observation space and an operator whose capability
  /// of simulating observations from this space is to be tested.
  oops::Parameter<std::vector<ObsTypeParameters<OBS>>> observations{"observations", {}, this};
};

// -----------------------------------------------------------------------------

/// \brief tests constructor and print method
template <typename OBS> void testConstructor() {
  typedef oops::ObsOperator<OBS>             ObsOperator_;
  typedef typename ObsOperator_::Parameters_ ObsOperatorParameters_;
  typedef ObsTypeParameters<OBS>             ObsTypeParameters_;
  typedef ObsTestsFixture<OBS>               Test_;
  typedef TestParameters<OBS>                TestParameters_;

  TestParameters_ testParams;
  testParams.validateAndDeserialize(TestEnvironment::config());

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const ObsTypeParameters_ &obsTypeParams = testParams.observations.value()[jj];
    const ObsOperatorParameters_ &obsOpParams = obsTypeParams.observer.obsOperator.value();
    if (obsTypeParams.expectConstructorToThrow.value() == boost::none) {
      auto hop = std::make_unique<ObsOperator_>(Test_::obspace()[jj], obsOpParams);;
      EXPECT(hop.get());
      oops::Log::info() << "Testing ObsOperator: " << *hop << std::endl;
      hop.reset();
      EXPECT(!hop.get());
    } else {
      // The constructor is expected to throw an exception containing the specified string.
      const std::string &expectedMessage =
          *obsTypeParams.expectConstructorToThrow.value();
      EXPECT_THROWS_MSG(ObsOperator_(Test_::obspace()[jj], obsOpParams),
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
  typedef ObsTypeParameters<OBS>       ObsTypeParameters_;
  typedef oops::ObsVector<OBS>         ObsVector_;
  typedef ObsTestsFixture<OBS>         Test_;
  typedef TestParameters<OBS>          TestParameters_;

  TestParameters_ testParams;
  testParams.validateAndDeserialize(TestEnvironment::config());

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const ObsTypeParameters_ &obsTypeParams = testParams.observations.value()[jj];

    if (obsTypeParams.expectConstructorToThrow.value() != boost::none)
      continue;

    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    ObsOperator_ hop(Test_::obspace()[jj], obsTypeParams.observer.obsOperator);

    // initialize bias correction
    const ObsAuxCtrl_ ybias(Test_::obspace()[jj], obsTypeParams.obsBias);

    // initialize geovals
    oops::Variables hopvars = hop.requiredVars();
    oops::Variables reducedHopvars = ybias.requiredVars();
    hopvars += reducedHopvars;  // the reduced format is derived from the sampled format
    // read geovals from the file (in the sampled format)
    GeoVaLs_ gval(obsTypeParams.geovals, Test_::obspace()[jj], hopvars);
    // convert geovals to the reduced format
    hop.computeReducedVars(reducedHopvars, gval);

    // create obsvector to hold H(x)
    ObsVector_ hofx(Test_::obspace()[jj]);

    // create obsvector to hold bias
    ObsVector_ bias(Test_::obspace()[jj]);
    bias.zero();

    // create diagnostics to hold HofX diags
    oops::Variables diagvars;
    diagvars += ybias.requiredHdiagnostics();
    ObsDiags_ diags(Test_::obspace()[jj], hop.locations(), diagvars);

    // call H(x), save result in the output file as @hofx
    if (obsTypeParams.expectSimulateObsToThrow.value() != boost::none) {
      // The simulateObs method is expected to throw an exception
      // containing the specified string.
      const std::string expectedMessage =
        *obsTypeParams.expectSimulateObsToThrow.value();
      EXPECT_THROWS_MSG(hop.simulateObs(gval, hofx, ybias, bias, diags),
                        expectedMessage.c_str());
      continue;
    } else {
      hop.simulateObs(gval, hofx, ybias, bias, diags);
    }
    hofx.save("hofx");
    bias.save("ObsBias");

    const double tol = obsTypeParams.tolerance.value().value();
    if (obsTypeParams.vectorRef.value() != boost::none) {
      // if reference h(x) is saved in file as a vector, read from file
      // and compare the norm of difference to zero
      ObsVector_ obsref(Test_::obspace()[jj], *obsTypeParams.vectorRef.value());
      obsref -= hofx;
      const double zz = obsref.rms();
      oops::Log::info() << "Vector difference between reference and computed: " << obsref;
      EXPECT(zz < 100*tol);  //  change tol from percent to actual value.
                             //  tol used in is_close is relative
    } else if (obsTypeParams.normRef.value() != boost::none) {
      // if reference h(x) is saved in file as a vector, read from file
      // and compare the difference, normalised by the reference values to zero
      ObsVector_ obsref(Test_::obspace()[jj], *obsTypeParams.normRef.value());
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
      const double xx = *obsTypeParams.rmsRef.value();

      oops::Log::debug() << "zz: " << std::fixed << std::setprecision(8) << zz << std::endl;
      oops::Log::debug() << "xx: " << std::fixed << std::setprecision(8) << xx << std::endl;

      EXPECT(oops::is_close(xx, zz, tol));
    }
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsOperator : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;

 public:
  ObsOperator() {
    // Needed because oops::ObsTypeParametersBase contains obs error parameters.
    oops::instantiateObsErrorFactory<OBS>();
  }

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
