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


#include "eckit/testing/Test.h"
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

const char *expectConstructorToThrow = "expect constructor to throw exception with message";
const char *expectSetTrajectoryToThrow = "expect setTrajectory to throw exception with message";
const char *expectSimulateObsToThrow = "expect simulateObs to throw exception with message";
const char *expectSimulateObsTLToThrow = "expect simulateObsTL to throw exception with message";
const char *expectSimulateObsADToThrow = "expect simulateObsAD to throw exception with message";

// -----------------------------------------------------------------------------
/// \brief tests constructor and print method
template <typename OBS> void testConstructor() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::LinearObsOperator<OBS>  LinearObsOperator_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::LocalConfiguration & conf = Test_::config(jj);
    // Use ObsOperator section of yaml unless LinearObsOperator is specified
    std::string confname = "obs operator";
    if (conf.has("linear obs operator")) confname = "linear obs operator";
    eckit::LocalConfiguration linobsopconf(conf, confname);

    if (!Test_::config(jj).has(expectConstructorToThrow)) {
      std::unique_ptr<LinearObsOperator_> linobsop(
        new LinearObsOperator_(Test_::obspace()[jj], linobsopconf));
      EXPECT(linobsop.get());
      oops::Log::test() << "Testing LinearObsOperator: " << *linobsop << std::endl;
      linobsop.reset();
      EXPECT(!linobsop.get());
    } else {
      // The constructor is expected to throw an exception containing the specified string.
      const std::string expectedMessage = Test_::config(jj).getString(expectConstructorToThrow);
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

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::LocalConfiguration & conf = Test_::config(jj);
    if (conf.has(expectConstructorToThrow) ||
        conf.has(expectSetTrajectoryToThrow) ||
        conf.has(expectSimulateObsToThrow) ||
        conf.has(expectSimulateObsTLToThrow))
      continue;

    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    eckit::LocalConfiguration obsopconf(conf, "obs operator");
    ObsOperator_ hop(Test_::obspace()[jj], obsopconf);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    // Use ObsOperator section of yaml unless LinearObsOperator is specified
    std::string confname = "obs operator";
    if (conf.has("linear obs operator")) confname = "linear obs operator";
    eckit::LocalConfiguration linobsopconf(conf, confname);
    LinearObsOperator_ hoptl(Test_::obspace()[jj], linobsopconf);

    // initialize obs bias
    eckit::LocalConfiguration biasconf = conf.getSubConfiguration("obs bias");
    typename ObsAuxCtrl_::Parameters_ biasparams;
    biasparams.validateAndDeserialize(biasconf);
    const ObsAuxCtrl_ ybias(Test_::obspace()[jj], biasparams);
    ObsAuxIncr_ ybinc(Test_::obspace()[jj], biasparams);

    // read geovals from the file
    const eckit::LocalConfiguration gconf(conf, "geovals");
    oops::Variables hopvars = hop.requiredVars();
    hopvars += ybias.requiredVars();
    const GeoVaLs_ gval(gconf, Test_::obspace()[jj], hopvars);

     // initialize Obs. Bias Covariance
    const ObsAuxCov_ Bobsbias(Test_::obspace()[jj], biasparams);

    // set trajectory for TL/AD to be the geovals from the file
    hoptl.setTrajectory(gval, ybias);

    // create obsvector
    ObsVector_ dy1(Test_::obspace()[jj]);

    // create geovals
    GeoVaLs_ dx(gconf, Test_::obspace()[jj], hoptl.requiredVars());

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

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::LocalConfiguration & conf = Test_::config(jj);
    if (conf.has(expectConstructorToThrow) ||
        conf.has(expectSetTrajectoryToThrow) ||
        conf.has(expectSimulateObsToThrow) ||
        conf.has(expectSimulateObsTLToThrow) ||
        conf.has(expectSimulateObsADToThrow))
      continue;

    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    eckit::LocalConfiguration obsopconf(conf, "obs operator");
    ObsOperator_ hop(Test_::obspace()[jj], obsopconf);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    // Use ObsOperator section of yaml unless LinearObsOperator is specified
    std::string confname = "obs operator";
    if (conf.has("linear obs operator")) confname = "linear obs operator";
    eckit::LocalConfiguration linobsopconf(conf, confname);
    LinearObsOperator_ hoptl(Test_::obspace()[jj], linobsopconf);

    const double tol = conf.getDouble("linear obs operator test.tolerance AD");
    // initialize bias correction
    eckit::LocalConfiguration biasconf = conf.getSubConfiguration("obs bias");
    typename ObsAuxCtrl_::Parameters_ biasparams;
    biasparams.validateAndDeserialize(biasconf);
    const ObsAuxCtrl_ ybias(Test_::obspace()[jj], biasparams);
    ObsAuxIncr_ ybinc1(Test_::obspace()[jj], biasparams);  // TL
    ObsAuxIncr_ ybinc2(Test_::obspace()[jj], biasparams);  // AD

    // initialize Obs. Bias Covariance
    const ObsAuxCov_ Bobsbias(Test_::obspace()[jj], biasparams);

    // read geovals from the file
    eckit::LocalConfiguration gconf(conf, "geovals");
    oops::Variables hopvars = hop.requiredVars();
    hopvars += ybias.requiredVars();
    const GeoVaLs_ gval(gconf, Test_::obspace()[jj], hopvars);

    // set TL/AD trajectory to the geovals from the file
    hoptl.setTrajectory(gval, ybias);

    ObsVector_ dy1(Test_::obspace()[jj]);
    ObsVector_ dy2(Test_::obspace()[jj]);
    GeoVaLs_ dx1(gconf, Test_::obspace()[jj], hoptl.requiredVars());
    GeoVaLs_ dx2(gconf, Test_::obspace()[jj], hoptl.requiredVars());

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

    oops::Log::test() << "Adjoint test result: (<x,HTy>-<Hx,y>)/<Hx,y> = "
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
  typedef oops::LinearObsOperator<OBS> LinearObsOperator_;
  typedef oops::ObsOperator<OBS>       ObsOperator_;
  typedef oops::ObsVector<OBS>         ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::LocalConfiguration & conf = Test_::config(jj);
    if (conf.has(expectConstructorToThrow) ||
        conf.has(expectSetTrajectoryToThrow) ||
        conf.has(expectSimulateObsToThrow) ||
        conf.has(expectSimulateObsTLToThrow))
      continue;

    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    eckit::LocalConfiguration obsopconf(conf, "obs operator");
    ObsOperator_ hop(Test_::obspace()[jj], obsopconf);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    // Use ObsOperator section of yaml unless LinearObsOperator is specified
    std::string confname = "obs operator";
    if (conf.has("linear obs operator")) confname = "linear obs operator";
    eckit::LocalConfiguration linobsopconf(conf, confname);
    LinearObsOperator_ hoptl(Test_::obspace()[jj], linobsopconf);

    const double tol = conf.getDouble("linear obs operator test.tolerance TL");
    const double alpha = conf.getDouble("linear obs operator test.coef TL", 0.1);
    const int iter = conf.getInt("linear obs operator test.iterations TL", 1);

    // initialize obs bias from file
    eckit::LocalConfiguration biasconf = conf.getSubConfiguration("obs bias");
    typename ObsAuxCtrl_::Parameters_ biasparams;
    biasparams.validateAndDeserialize(biasconf);
    const ObsAuxCtrl_ ybias0(Test_::obspace()[jj], biasparams);
    ObsAuxCtrl_ ybias(Test_::obspace()[jj], biasparams);

    // initialize Obs. Bias Covariance
    const ObsAuxCov_ Bobsbias(Test_::obspace()[jj], biasparams);

    // read geovals from the file
    const eckit::LocalConfiguration gconf(conf, "geovals");
    oops::Variables hopvars = hop.requiredVars();
    hopvars += ybias0.requiredVars();
    const GeoVaLs_ x0(gconf, Test_::obspace()[jj], hopvars);
    GeoVaLs_ x(gconf, Test_::obspace()[jj], hopvars);

    // set TL trajectory to the geovals and the bias coeff. from the files
    hoptl.setTrajectory(x0, ybias0);

    // create obsvectors
    ObsVector_ y1(Test_::obspace()[jj]);
    ObsVector_ y2(Test_::obspace()[jj]);
    ObsVector_ y3(Test_::obspace()[jj]);

    // create obsdatavector to hold diags
    oops::Variables diagvars;
    diagvars += ybias0.requiredHdiagnostics();
    ObsDiags_ ydiag(Test_::obspace()[jj], hop.locations(), diagvars);

    // y1 = hop(x0, ybias0)
    hop.simulateObs(x0, y1, ybias0, ydiag);

    // randomize dx and ybinc
    GeoVaLs_ dx(gconf, Test_::obspace()[jj], hoptl.requiredVars());
    dx.random();
    ObsAuxIncr_ ybinc(Test_::obspace()[jj], biasparams);
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

      // y2 = hop(x0+alpha*dx, ybias0+alpha*ybinc)
      hop.simulateObs(x, y2, ybias, ydiag);
      y2 -= y1;
      // y3 = hoptl(alpha*dx, alpha*ybinc)
      hoptl.simulateObsTL(dx, y3, ybinc);
      y2 -= y3;

      double test_norm = y2.rms();
      oops::Log::test() << "Iter:" << jter << " ||(h(x+alpha*dx)-h(x)-h'*(alpha*dx))||="
                        << test_norm << std::endl;
    }
    EXPECT(y2.rms() < tol);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testException() {
  typedef ObsTestsFixture<OBS> Test_;
  typedef oops::GeoVaLs<OBS>           GeoVaLs_;
  typedef oops::ObsOperator<OBS>       ObsOperator_;
  typedef oops::LinearObsOperator<OBS> LinearObsOperator_;
  typedef oops::ObsAuxControl<OBS>     ObsAuxCtrl_;
  typedef oops::ObsAuxIncrement<OBS>   ObsAuxIncr_;
  typedef oops::ObsAuxCovariance<OBS>  ObsAuxCov_;
  typedef oops::ObsVector<OBS>         ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const eckit::LocalConfiguration & conf = Test_::config(jj);
    if (conf.has(expectConstructorToThrow))
      continue;

    // Set up objects prior to throwing exceptions.
    eckit::LocalConfiguration obsopconf(conf, "obs operator");
    ObsOperator_ hop(Test_::obspace()[jj], obsopconf);
    std::string confname = "obs operator";
    if (conf.has("linear obs operator")) confname = "linear obs operator";
    eckit::LocalConfiguration linobsopconf(conf, confname);
    LinearObsOperator_ hoptl(Test_::obspace()[jj], linobsopconf);
    eckit::LocalConfiguration biasconf = conf.getSubConfiguration("obs bias");
    typename ObsAuxCtrl_::Parameters_ biasparams;
    biasparams.validateAndDeserialize(biasconf);
    const ObsAuxCtrl_ ybias(Test_::obspace()[jj], biasparams);
    ObsAuxIncr_ ybinc(Test_::obspace()[jj], biasparams);
    const ObsAuxCov_ Bobsbias(Test_::obspace()[jj], biasparams);
    eckit::LocalConfiguration gconf(conf, "geovals");
    oops::Variables hopvars = hop.requiredVars();
    hopvars += ybias.requiredVars();
    const GeoVaLs_ gval(gconf, Test_::obspace()[jj], hopvars);
    oops::Variables diagvars;
    diagvars += ybias.requiredHdiagnostics();

    if (Test_::config(jj).has(expectSetTrajectoryToThrow)) {
      // The setTrajectory method is expected to throw an exception
      // containing the specified string.
      const std::string expectedMessage =
        Test_::config(jj).getString(expectSetTrajectoryToThrow);
      EXPECT_THROWS_MSG(hoptl.setTrajectory(gval, ybias),
                        expectedMessage.c_str());
      // Do not continue further because setTrajectory must be run
      // before simulateObsTL and simulateObsAD.
      continue;
    }

    if (Test_::config(jj).has(expectSimulateObsTLToThrow)) {
      hoptl.setTrajectory(gval, ybias);
      ObsVector_ dy1(Test_::obspace()[jj]);
      GeoVaLs_ dx1(gconf, Test_::obspace()[jj], hoptl.requiredVars());
      dx1.random();
      Bobsbias.randomize(ybinc);
      // The simulateObsTL method is expected to throw an exception
      // containing the specified string.
      const std::string expectedMessage =
        Test_::config(jj).getString(expectSimulateObsTLToThrow);
      EXPECT_THROWS_MSG(hoptl.simulateObsTL(dx1, dy1, ybinc),
                        expectedMessage.c_str());
    }

    if (Test_::config(jj).has(expectSimulateObsADToThrow)) {
      hoptl.setTrajectory(gval, ybias);
      ObsVector_ dy2(Test_::obspace()[jj]);
      GeoVaLs_ dx2(gconf, Test_::obspace()[jj], hoptl.requiredVars());
      Bobsbias.randomize(ybinc);
      dy2.random();
      dx2.zero();
      ybinc.zero();
      // The simulateObsAD method is expected to throw an exception
      // containing the specified string.
      const std::string expectedMessage =
        Test_::config(jj).getString(expectSimulateObsADToThrow);
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
  LinearObsOperator() {}
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
