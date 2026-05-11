/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_LINEAROBSOPERATOR_H_
#define TEST_INTERFACE_LINEAROBSOPERATOR_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0


#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxCovariance.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/runs/Test.h"
#include "oops/util/dot_product.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

/// \brief Extract linear obs operator config from the 'linear obs operator' YAML option
/// if present or from the 'obs operator' option otherwise.
eckit::LocalConfiguration linearObsOperatorConf(const eckit::Configuration & oconf) {
  if (oconf.has("linear obs operator"))
    return eckit::LocalConfiguration(oconf, "linear obs operator");
  else
    // [Comment copied from ObserverTLAD.h]
    // Hack: when "linear obs operator" is not specified in the input file, reinterpret
    //       the entry for "obs operator" as a linear obs operator option. In the long
    //       term, we need a design that either,
    //       - allows constructing LinearObsOperator from either set of Parameters, or
    //       - merges the two sets of Parameters so this switch can be removed
    return eckit::LocalConfiguration(oconf, "obs operator");
}

// -----------------------------------------------------------------------------
/// \brief tests constructor and print method
template <typename OBS> void testConstructor() {
  typedef oops::LinearObsOperator<OBS>  LinearObsOperator_;
  typedef ObsTestsFixture<OBS>          Test_;

  const std::vector<eckit::LocalConfiguration> obsConfs =
      TestEnvironment::config().getSubConfigurations("observations");

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::LocalConfiguration linobsopconf = linearObsOperatorConf(obsConfs[jj]);

    if (!obsConfs[jj].has("expect constructor to throw exception with message")) {
      std::unique_ptr<LinearObsOperator_> linobsop(
        new LinearObsOperator_(Test_::obspace()[jj], linobsopconf));
      EXPECT(linobsop.get());
      oops::Log::info() << "Testing LinearObsOperator: " << *linobsop << std::endl;
      linobsop.reset();
      EXPECT(!linobsop.get());
    } else {
      // The constructor is expected to throw an exception containing the specified string.
      const std::string expectedMessage =
          obsConfs[jj].getString("expect constructor to throw exception with message");
      EXPECT_THROWS_MSG(LinearObsOperator_(Test_::obspace()[jj], linobsopconf),
                        expectedMessage.c_str());
    }
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testLinearity() {
  typedef ObsTestsFixture<OBS>         Test_;
  typedef oops::GeoVaLs<OBS>           GeoVaLs_;
  typedef oops::ObsAuxControl<OBS>     ObsAuxCtrl_;
  typedef oops::ObsAuxIncrement<OBS>   ObsAuxIncr_;
  typedef oops::ObsAuxCovariance<OBS>  ObsAuxCov_;
  typedef oops::ObsOperator<OBS>       ObsOperator_;
  typedef oops::LinearObsOperator<OBS> LinearObsOperator_;
  typedef oops::ObsVector<OBS>         ObsVector_;

  const double zero = 0.0;
  const double coef = 3.14;
  const double tol = 1.0e-11;

  const std::vector<eckit::LocalConfiguration> obsConfs =
      TestEnvironment::config().getSubConfigurations("observations");

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::Configuration & obsconf = obsConfs[jj];
    if (obsconf.has("expect constructor to throw exception with message") ||
        obsconf.has("expect setTrajectory to throw exception with message") ||
        obsconf.has("expect simulateObs to throw exception with message") ||
        obsconf.has("expect simulateObsTL to throw exception with message"))
      continue;

    const eckit::LocalConfiguration oopconf(obsconf, "obs operator");
    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    ObsOperator_ hop(Test_::obspace()[jj], oopconf);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    LinearObsOperator_ hoptl(Test_::obspace()[jj], linearObsOperatorConf(obsconf));

    // initialize qc flags ObsDataVector
    oops::ObsDataVector<OBS, int> qc_flags(
      Test_::obspace()[jj],
      Test_::obspace()[jj].obsvariables(),
      std::string());

    // Check if linear obs operator test config contains QCFlagsGroupName option
    // Read the group_name and from obs space read the values
    if (obsconf.has("linear obs operator test")) {
      const eckit::LocalConfiguration lotConf(obsconf, "linear obs operator test");
      if (lotConf.has("QCFlagsGroupName")) qc_flags.read(lotConf.getString("QCFlagsGroupName"));
    }

    // initialize obs bias
    const eckit::LocalConfiguration bconf = obsconf.getSubConfiguration("obs bias");
    const ObsAuxCtrl_ ybias(Test_::obspace()[jj], bconf);
    ObsAuxIncr_ ybinc(Test_::obspace()[jj], bconf);

    // initialize geovals
    oops::Variables hopvars = hop.requiredVars();
    oops::Variables reducedHopvars = ybias.requiredVars();
    hopvars += reducedHopvars;
    // read geovals from the file (in the sampled format)
    GeoVaLs_ gval(eckit::LocalConfiguration(obsconf, "geovals"), Test_::obspace()[jj], hopvars);
    // convert geovals to the reduced format
    hop.computeReducedVars(reducedHopvars, gval);

     // initialize Obs. Bias Covariance
    const ObsAuxCov_ Bobsbias(Test_::obspace()[jj], bconf);

    // set trajectory for TL/AD to be the geovals from the file
    hoptl.setTrajectory(gval, ybias, qc_flags);

    // create obsvector
    ObsVector_ dy1(Test_::obspace()[jj]);

    // create geovals
    const oops::Variables hoptlvars = hoptl.requiredVars();
    GeoVaLs_ dx(eckit::LocalConfiguration(obsconf, "geovals"), Test_::obspace()[jj], hoptlvars);

    // test rms(H * (dx, ybinc)) = 0, when dx = 0
    dx.zero();
    ybinc.zero();
    hoptl.simulateObsTL(dx, dy1, ybinc);
    EXPECT(dy1.rms() == zero);

    // test rms(H * (dx, ybinc)) > 0, when dx is random
    dx.random();
    Bobsbias.randomize(ybinc);
    hoptl.simulateObsTL(dx, dy1, ybinc);
    EXPECT(dy1.rms() > zero);

    // test k * H * (dx, ybinc) ~ H * (k*dx, k*ybinc)
    dy1 *= coef;
    dx  *= coef;
    ybinc *= coef;
    ObsVector_ dy2(Test_::obspace()[jj]);
    hoptl.simulateObsTL(dx, dy2, ybinc);

    dy2 -= dy1;
    EXPECT(dy2.rms() / dy1.rms() < tol);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testAdjoint() {
  typedef ObsTestsFixture<OBS> Test_;
  typedef oops::GeoVaLs<OBS>           GeoVaLs_;
  typedef oops::ObsOperator<OBS>       ObsOperator_;
  typedef oops::LinearObsOperator<OBS> LinearObsOperator_;
  typedef oops::ObsAuxControl<OBS>     ObsAuxCtrl_;
  typedef oops::ObsAuxIncrement<OBS>   ObsAuxIncr_;
  typedef oops::ObsAuxCovariance<OBS>  ObsAuxCov_;
  typedef oops::ObsVector<OBS>         ObsVector_;
  const double zero = 0.0;

  const std::vector<eckit::LocalConfiguration> obsConfs =
      TestEnvironment::config().getSubConfigurations("observations");

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::Configuration & obsconf = obsConfs[jj];
    if (obsconf.has("expect constructor to throw exception with message") ||
        obsconf.has("expect setTrajectory to throw exception with message") ||
        obsconf.has("expect simulateObs to throw exception with message") ||
        obsconf.has("expect simulateObsTL to throw exception with message") ||
        obsconf.has("expect simulateObsAD to throw exception with message"))
      continue;

    const eckit::LocalConfiguration oopconf(obsconf, "obs operator");
    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    ObsOperator_ hop(Test_::obspace()[jj], oopconf);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    LinearObsOperator_ hoptl(Test_::obspace()[jj], linearObsOperatorConf(obsconf));

    // initialize qc flags ObsDataVector
    oops::ObsDataVector<OBS, int> qc_flags(Test_::obspace()[jj],
                                           Test_::obspace()[jj].obsvariables(), std::string());

    // Check if linear obs operator test config contains QCFlagsGroupName option
    // Read the group_name and from obs space read the values
    const eckit::LocalConfiguration lotConf(obsconf, "linear obs operator test");
    if (lotConf.has("QCFlagsGroupName"))
      qc_flags.read(lotConf.getString("QCFlagsGroupName"));

    const double tol = lotConf.getDouble("tolerance AD");
    // initialize bias correction
    const eckit::LocalConfiguration bconf = obsconf.getSubConfiguration("obs bias");
    const ObsAuxCtrl_ ybias(Test_::obspace()[jj], bconf);
    ObsAuxIncr_ ybinc1(Test_::obspace()[jj], bconf);
    ObsAuxIncr_ ybinc2(Test_::obspace()[jj], bconf);

    // initialize Obs. Bias Covariance
    const ObsAuxCov_ Bobsbias(Test_::obspace()[jj], bconf);

    // initialize geovals
    oops::Variables hopvars = hop.requiredVars();
    oops::Variables reducedHopvars = ybias.requiredVars();
    hopvars += reducedHopvars;  // the reduced format is derived from the sampled format
    // read geovals from the file (in the sampled format)
    GeoVaLs_ gval(eckit::LocalConfiguration(obsconf, "geovals"), Test_::obspace()[jj], hopvars);
    // convert geovals to the reduced format
    hop.computeReducedVars(reducedHopvars, gval);

    // set TL/AD trajectory to the geovals from the file
    hoptl.setTrajectory(gval, ybias, qc_flags);

    ObsVector_ dy1(Test_::obspace()[jj]);
    ObsVector_ dy2(Test_::obspace()[jj]);
    const oops::Variables hoptlvars = hoptl.requiredVars();
    GeoVaLs_ dx1(eckit::LocalConfiguration(obsconf, "geovals"), Test_::obspace()[jj], hoptlvars);
    GeoVaLs_ dx2(eckit::LocalConfiguration(obsconf, "geovals"), Test_::obspace()[jj], hoptlvars);

    // calculate dy1 = H (dx1, ybinc1) (with random dx1, and random ybinc1)
    dx1.random();
    EXPECT(dot_product(dx1, dx1) > zero);  //  BOOST_REQUIRE
    Bobsbias.randomize(ybinc1);
    hoptl.simulateObsTL(dx1, dy1, ybinc1);
    EXPECT(dot_product(dy1, dy1) > zero);

    // calculate (dx2, ybinc2) = HT dy2 (with random dy2)
    dy2.random();
    EXPECT(dot_product(dy2, dy2) > zero);  //  BOOST_REQUIRE
    dx2.zero();
    ybinc2.zero();
    hoptl.simulateObsAD(dx2, dy2, ybinc2);
    EXPECT(dot_product(dx2, dx2) > zero);

    const double zz1 = dot_product(dx1, dx2) + dot_product(ybinc1, ybinc2);
    const double zz2 = dot_product(dy1, dy2);

    oops::Log::info() << "Adjoint test result: (<x,HTy>-<Hx,y>)/<Hx,y> = "
                       << (zz1-zz2)/zz2 << std::endl;

    EXPECT(zz1 != zero);
    EXPECT(zz2 != zero);
    EXPECT(oops::is_close(zz1, zz2, tol));
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testTangentLinear() {
  // Test  ||(hop(x+alpha*dx)-hop(x)) - hoptl(alpha*dx)|| < tol
  typedef ObsTestsFixture<OBS>         Test_;
  typedef oops::GeoVaLs<OBS>           GeoVaLs_;
  typedef oops::ObsDiagnostics<OBS>    ObsDiags_;
  typedef oops::ObsAuxControl<OBS>     ObsAuxCtrl_;
  typedef oops::ObsAuxIncrement<OBS>   ObsAuxIncr_;
  typedef oops::ObsAuxCovariance<OBS>  ObsAuxCov_;
  typedef oops::ObsOperator<OBS>       ObsOperator_;
  typedef oops::LinearObsOperator<OBS> LinearObsOperator_;
  typedef oops::ObsVector<OBS>         ObsVector_;

  const std::vector<eckit::LocalConfiguration> obsConfs =
      TestEnvironment::config().getSubConfigurations("observations");

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::Configuration & obsconf = obsConfs[jj];
    if (obsconf.has("expect constructor to throw exception with message") ||
        obsconf.has("expect setTrajectory to throw exception with message") ||
        obsconf.has("expect simulateObs to throw exception with message") ||
        obsconf.has("expect simulateObsTL to throw exception with message"))
      continue;

    const eckit::LocalConfiguration oopconf(obsconf, "obs operator");
    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    ObsOperator_ hop(Test_::obspace()[jj], oopconf);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    LinearObsOperator_ hoptl(Test_::obspace()[jj], linearObsOperatorConf(obsconf));

    const eckit::LocalConfiguration lotConf(obsconf, "linear obs operator test");
    const double tol = lotConf.getDouble("tolerance TL");
    const double alpha = lotConf.getDouble("coef TL", 0.1);
    const int iter = lotConf.getInt("iterations TL", 1);

    // initialize obs bias from file
    const eckit::LocalConfiguration bconf = obsconf.getSubConfiguration("obs bias");
    const ObsAuxCtrl_ ybias0(Test_::obspace()[jj], bconf);
    ObsAuxCtrl_ ybias(Test_::obspace()[jj], bconf);

    // initialize Obs. Bias Covariance
    const ObsAuxCov_ Bobsbias(Test_::obspace()[jj], bconf);

    // initialize geovals
    oops::Variables hopvars = hop.requiredVars();
    oops::Variables reducedHopvars = ybias0.requiredVars();
    hopvars += reducedHopvars;  // the reduced format is derived from the sampled format
    // read geovals from the file
    const eckit::LocalConfiguration geovalsConf(obsconf, "geovals");
    GeoVaLs_ x0(geovalsConf, Test_::obspace()[jj], hopvars);
    GeoVaLs_ x(geovalsConf, Test_::obspace()[jj], hopvars);
    // convert geovals to the reduced format
    hop.computeReducedVars(reducedHopvars, x0);
    hop.computeReducedVars(reducedHopvars, x);

    // create obsvectors
    ObsVector_ y1(Test_::obspace()[jj]);
    ObsVector_ y2(Test_::obspace()[jj]);
    ObsVector_ y3(Test_::obspace()[jj]);
    ObsVector_ bias(Test_::obspace()[jj]);

    // initialize qc flags ObsDataVector
    oops::ObsDataVector<OBS, int> qc_flags(Test_::obspace()[jj],
                                           Test_::obspace()[jj].obsvariables(), std::string());

    // Check if linear obs operator test config contains QCFlagsGroupName option
    // Read the group_name and from obs space read the values
    if (lotConf.has("QCFlagsGroupName")) qc_flags.read(lotConf.getString("QCFlagsGroupName"));

    // set TL trajectory to the geovals and the bias coeff. from the files
    hoptl.setTrajectory(x0, ybias0, qc_flags);

    bias.zero();

    // create obsdatavector to hold diags
    oops::ObsVariables diagvars;
    diagvars += ybias0.requiredHdiagnostics();
    ObsDiags_ ydiag(Test_::obspace()[jj], hop.locations(), diagvars);

    // y1 = hop(x0, ybias0)
    hop.simulateObs(x0, y1, ybias0, qc_flags, bias, ydiag);

    // randomize dx and ybinc
    const oops::Variables hoptlvars = hoptl.requiredVars();
    GeoVaLs_ dx(geovalsConf, Test_::obspace()[jj], hoptlvars);
    dx.random();
    ObsAuxIncr_ ybinc(Test_::obspace()[jj], bconf);
    Bobsbias.randomize(ybinc);

    // scale dx by x0
    dx *= x0;

    for (int jter = 0; jter < iter; ++jter) {
      // x = x0 + alpha*dx
      dx *= alpha;
      x = x0;
      x += dx;
      // ybias = ybias0 + alpha*ybinc
      ybinc *= alpha;
      ybias = ybias0;
      ybias += ybinc;
      bias.zero();

      // y2 = hop(x0+alpha*dx, ybias0+alpha*ybinc)
      hop.simulateObs(x, y2, ybias, qc_flags, bias, ydiag);
      y2 -= y1;
      // y3 = hoptl(alpha*dx, alpha*ybinc)
      hoptl.simulateObsTL(dx, y3, ybinc);
      y2 -= y3;

      double test_norm = y2.rms();
      oops::Log::info() << "Iter:" << jter << " ||(h(x+alpha*dx)-h(x)-h'*(alpha*dx))||="
                        << test_norm << std::endl;
    }
    EXPECT(y2.rms() < tol);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testException() {
  typedef ObsTestsFixture<OBS>         Test_;
  typedef oops::GeoVaLs<OBS>           GeoVaLs_;
  typedef oops::ObsOperator<OBS>       ObsOperator_;
  typedef oops::LinearObsOperator<OBS> LinearObsOperator_;
  typedef oops::ObsAuxControl<OBS>     ObsAuxCtrl_;
  typedef oops::ObsAuxIncrement<OBS>   ObsAuxIncr_;
  typedef oops::ObsAuxCovariance<OBS>  ObsAuxCov_;
  typedef oops::ObsVector<OBS>         ObsVector_;

  const std::vector<eckit::LocalConfiguration> obsConfs =
      TestEnvironment::config().getSubConfigurations("observations");

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::Configuration & obsconf = obsConfs[jj];
    if (obsconf.has("expect constructor to throw exception with message"))
      continue;

    // Set up objects prior to throwing exceptions.
    const eckit::LocalConfiguration oopconf(obsconf, "obs operator");
    ObsOperator_ hop(Test_::obspace()[jj], oopconf);
    LinearObsOperator_ hoptl(Test_::obspace()[jj], linearObsOperatorConf(obsconf));
    const eckit::LocalConfiguration bconf = obsconf.getSubConfiguration("obs bias");
    const ObsAuxCtrl_ ybias(Test_::obspace()[jj], bconf);
    ObsAuxIncr_ ybinc(Test_::obspace()[jj], bconf);
    const ObsAuxCov_ Bobsbias(Test_::obspace()[jj], bconf);
    oops::Variables hopvars = hop.requiredVars();
    oops::Variables reducedHopvars = ybias.requiredVars();
    hopvars += reducedHopvars;
    GeoVaLs_ gval(eckit::LocalConfiguration(obsconf, "geovals"), Test_::obspace()[jj], hopvars);
    hop.computeReducedVars(reducedHopvars, gval);
    oops::ObsVariables diagvars;
    diagvars += ybias.requiredHdiagnostics();
    const oops::Variables hoptlvars = hoptl.requiredVars();

    // initialize qc flags ObsDataVector
    oops::ObsDataVector<OBS, int> qc_flags(Test_::obspace()[jj],
                                           Test_::obspace()[jj].obsvariables(), std::string());

    // Check if linear obs operator test config contains QCFlagsGroupName option
    // Read the group_name and from obs space read the values
    if (obsconf.has("linear obs operator test")) {
      const eckit::LocalConfiguration lotConf(obsconf, "linear obs operator test");
      if (lotConf.has("QCFlagsGroupName")) qc_flags.read(lotConf.getString("QCFlagsGroupName"));
    }

    if (obsconf.has("expect setTrajectory to throw exception with message")) {
      // The setTrajectory method is expected to throw an exception
      // containing the specified string.
      const std::string expectedMessage =
          obsconf.getString("expect setTrajectory to throw exception with message");
      EXPECT_THROWS_MSG(hoptl.setTrajectory(gval, ybias, qc_flags),
                        expectedMessage.c_str());
      // Do not continue further because setTrajectory must be run
      // before simulateObsTL and simulateObsAD.
      continue;
    }
    if (obsconf.has("expect simulateObsTL to throw exception with message")) {
      hoptl.setTrajectory(gval, ybias, qc_flags);
      ObsVector_ dy1(Test_::obspace()[jj]);
      GeoVaLs_ dx1(eckit::LocalConfiguration(obsconf, "geovals"), Test_::obspace()[jj], hoptlvars);
      dx1.random();
      Bobsbias.randomize(ybinc);
      // The simulateObsTL method is expected to throw an exception
      // containing the specified string.
      const std::string expectedMessage =
          obsconf.getString("expect simulateObsTL to throw exception with message");
      EXPECT_THROWS_MSG(hoptl.simulateObsTL(dx1, dy1, ybinc),
                        expectedMessage.c_str());
    }

    if (obsconf.has("expect simulateObsAD to throw exception with message")) {
      hoptl.setTrajectory(gval, ybias, qc_flags);
      ObsVector_ dy2(Test_::obspace()[jj]);
      GeoVaLs_ dx2(eckit::LocalConfiguration(obsconf, "geovals"), Test_::obspace()[jj], hoptlvars);
      Bobsbias.randomize(ybinc);
      dy2.random();
      dx2.zero();
      ybinc.zero();
      // The simulateObsAD method is expected to throw an exception
      // containing the specified string.
      const std::string expectedMessage =
          obsconf.getString("expect simulateObsAD to throw exception with message");
      EXPECT_THROWS_MSG(hoptl.simulateObsAD(dx2, dy2, ybinc),
                        expectedMessage.c_str());
    }
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
class LinearObsOperator : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;

 public:
  explicit LinearObsOperator(const eckit::mpi::Comm & comm = oops::mpi::world()) :
    oops::Test(comm) {}
  virtual ~LinearObsOperator() {}

 private:
  std::string testid() const override {return "test::LinearObsOperator<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/LinearObsOperator/testConstructor")
      { testConstructor<OBS>(); });
    ts.emplace_back(CASE("interface/LinearObsOperator/testLinearity")
      { testLinearity<OBS>(); });
    ts.emplace_back(CASE("interface/LinearObsOperator/testTangentLinear")
      { testTangentLinear<OBS>(); });
    ts.emplace_back(CASE("interface/LinearObsOperator/testAdjoint")
      { testAdjoint<OBS>(); });
    ts.emplace_back(CASE("interface/LinearObsOperator/testException")
      { testException<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LINEAROBSOPERATOR_H_
