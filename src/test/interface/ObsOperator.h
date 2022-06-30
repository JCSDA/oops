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

const char *expectConstructorToThrow = "expect constructor to throw exception with message";
const char *expectSimulateObsToThrow = "expect simulateObs to throw exception with message";

// -----------------------------------------------------------------------------
/// \brief tests constructor and print method
template <typename OBS> void testConstructor() {
  typedef ObsTestsFixture<OBS> Test_;
  typedef oops::ObsOperator<OBS>       ObsOperator_;
  typedef typename ObsOperator_::Parameters_ ObsOperatorParameters_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration obsopconf(Test_::config(jj), "obs operator");
    ObsOperatorParameters_ obsopparams;
    obsopparams.validateAndDeserialize(obsopconf);

    if (!Test_::config(jj).has(expectConstructorToThrow)) {
      std::unique_ptr<ObsOperator_> hop(new ObsOperator_(Test_::obspace()[jj], obsopparams));
      EXPECT(hop.get());
      oops::Log::test() << "Testing ObsOperator: " << *hop << std::endl;
      hop.reset();
      EXPECT(!hop.get());
    } else {
      // The constructor is expected to throw an exception containing the specified string.
      const std::string expectedMessage = Test_::config(jj).getString(expectConstructorToThrow);
      EXPECT_THROWS_MSG(ObsOperator_(Test_::obspace()[jj], obsopparams),
                        expectedMessage.c_str());
    }
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testSimulateObs() {
  typedef ObsTestsFixture<OBS> Test_;
  typedef oops::GeoVaLs<OBS>           GeoVaLs_;
  typedef typename GeoVaLs_::Parameters_ GeoVaLsParameters_;
  typedef oops::ObsDiagnostics<OBS>    ObsDiags_;
  typedef oops::ObsAuxControl<OBS>     ObsAuxCtrl_;
  typedef oops::ObsOperator<OBS>       ObsOperator_;
  typedef typename ObsOperator_::Parameters_ ObsOperatorParameters_;
  typedef oops::ObsVector<OBS>         ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::LocalConfiguration & conf = Test_::config(jj);
    if (conf.has(expectConstructorToThrow))
      continue;

    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    eckit::LocalConfiguration obsopconf(conf, "obs operator");
    ObsOperatorParameters_ obsopparams;
    obsopparams.validateAndDeserialize(obsopconf);
    ObsOperator_ hop(Test_::obspace()[jj], obsopparams);

    // initialize bias correction
    eckit::LocalConfiguration biasconf = conf.getSubConfiguration("obs bias");
    typename ObsAuxCtrl_::Parameters_ biasparams;
    biasparams.validateAndDeserialize(biasconf);
    const ObsAuxCtrl_ ybias(Test_::obspace()[jj], biasparams);

    // read geovals from the file
    eckit::LocalConfiguration gconf(conf, "geovals");
    GeoVaLsParameters_ geovalsparams;
    geovalsparams.validateAndDeserialize(gconf);
    oops::Variables hopvars = hop.requiredVars();
    hopvars += ybias.requiredVars();
    const GeoVaLs_ gval(geovalsparams, Test_::obspace()[jj], hopvars);

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
    if (Test_::config(jj).has(expectSimulateObsToThrow)) {
      // The simulateObs method is expected to throw an exception
      // containing the specified string.
      const std::string expectedMessage =
        Test_::config(jj).getString(expectSimulateObsToThrow);
      EXPECT_THROWS_MSG(hop.simulateObs(gval, hofx, ybias, bias, diags),
                        expectedMessage.c_str());
      continue;
    } else {
      hop.simulateObs(gval, hofx, ybias, bias, diags);
    }
    hofx.save("hofx");
    bias.save("ObsBias");

    const double tol = conf.getDouble("tolerance");
    if (conf.has("vector ref")) {
      // if reference h(x) is saved in file as a vector, read from file
      // and compare the norm of difference to zero
      ObsVector_ obsref(Test_::obspace()[jj], conf.getString("vector ref"));
      obsref -= hofx;
      const double zz = obsref.rms();
      oops::Log::test() << "Vector difference between reference and computed: " << obsref;
      EXPECT(zz < 100*tol);  //  change tol from percent to actual value.
                             //  tol used in is_close is relative
    } else if (conf.has("norm ref")) {
      // if reference h(x) is saved in file as a vector, read from file
      // and compare the difference, normalised by the reference values to zero
      ObsVector_ obsref(Test_::obspace()[jj], conf.getString("norm ref"));
      obsref -= hofx;
      obsref /= hofx;
      const double zz = obsref.rms();
      oops::Log::test() << "Normalised vector difference between reference and computed: "
                        << obsref;
      EXPECT(zz < 100*tol);  //  change tol from percent to actual value.
                             //  tol used in is_close is relative
    } else {
      // else compare h(x) norm to the norm from the config
      const double zz = hofx.rms();
      const double xx = conf.getDouble("rms ref");

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
  ObsOperator() {}
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
