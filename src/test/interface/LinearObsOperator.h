/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_LINEAROBSOPERATOR_H_
#define TEST_INTERFACE_LINEAROBSOPERATOR_H_

#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/scoped_ptr.hpp>

#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsOperator.h"
#include "oops/runs/Test.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ObsTestsFixture<MODEL>  Test_;
  typedef oops::LinearObsOperator<MODEL>  LinearObsOperator_;

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    boost::scoped_ptr<LinearObsOperator_> ov(
      new LinearObsOperator_(Test_::obspace()[jj], conf[jj]));
    BOOST_CHECK(ov.get());

    ov.reset();
    BOOST_CHECK(!ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearity() {
  typedef ObsTestsFixture<MODEL>         Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::ObsAuxIncrement<MODEL>   ObsAuxIncr_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;
  typedef oops::LinearObsOperator<MODEL> LinearObsOperator_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

  const double zero = 0.0;
  const double coef = 3.14;
  const double tol = TestEnvironment::config().getDouble("LinearObsOpTest.toleranceAD");
  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    ObsOperator_ hop(Test_::obspace()[jj], conf[jj]);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    LinearObsOperator_ hoptl(Test_::obspace()[jj], conf[jj]);

    // read geovals from the file
    const eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    const GeoVaLs_ gval(gconf, hop.variables());

    // initialize obs bias
    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);
    const ObsAuxIncr_ ybinc(biasConf);

    // set trajectory for TL/AD to be the geovals from the file
    hoptl.setTrajectory(gval, ybias);

    // create obsvector
    ObsVector_ dy1(Test_::obspace()[jj], hop.observed());

    // create geovals
    GeoVaLs_ dx(gconf, hoptl.variables());

    // test rms(Hdx) = 0, when dx = 0
    dx.zero();
    hoptl.simulateObsTL(dx, dy1, ybinc);
    BOOST_CHECK_EQUAL(dy1.rms(), zero);

    // test rms(Hdx) > 0, when dx is random
    dx.random();
    hoptl.simulateObsTL(dx, dy1, ybinc);
    BOOST_CHECK(dy1.rms() > zero);

    // test k * H * dx ~ H * (k*dx)
    dy1 *= coef;
    dx  *= coef;
    ObsVector_ dy2(Test_::obspace()[jj], hop.observed());
    hoptl.simulateObsTL(dx, dy2, ybinc);

    dy1 -= dy2;

    BOOST_CHECK_SMALL(dy1.rms(), tol);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testAdjoint() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;
  typedef oops::LinearObsOperator<MODEL> LinearObsOperator_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::ObsAuxIncrement<MODEL>   ObsAuxIncr_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

  const double zero = 0.0;
  const double tol = TestEnvironment::config().getDouble("LinearObsOpTest.toleranceAD");
  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    ObsOperator_ hop(Test_::obspace()[jj], conf[jj]);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    LinearObsOperator_ hoptl(Test_::obspace()[jj], conf[jj]);

    // read geovals from the file
    eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    const GeoVaLs_ gval(gconf, hop.variables());

    // initialize bias correction
    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);
    ObsAuxIncr_ ybinc(biasConf);

    // set TL/AD trajectory to the geovals from the file
    hoptl.setTrajectory(gval, ybias);

    ObsVector_ dy1(Test_::obspace()[jj], hop.observed());
    ObsVector_ dy2(Test_::obspace()[jj], hop.observed());
    GeoVaLs_ dx1(gconf, hoptl.variables());
    GeoVaLs_ dx2(gconf, hoptl.variables());

    // calculate dy1 = H dx1 (with random dx1)
    dx1.random();
    BOOST_REQUIRE(dot_product(dx1, dx1) > zero);
    hoptl.simulateObsTL(dx1, dy1, ybinc);
    BOOST_CHECK(dot_product(dy1, dy1) > zero);

    // calculate dx2 = HT dy2 (with random dy2)
    dy2.random();
    BOOST_REQUIRE(dot_product(dy2, dy2) > zero);
    dx2.zero();
    hoptl.simulateObsAD(dx2, dy2, ybinc);
    BOOST_CHECK(dot_product(dx2, dx2) > zero);

    const double zz1 = dot_product(dx1, dx2);
    const double zz2 = dot_product(dy1, dy2);

    oops::Log::info() << "Adjoint test result: (<x,HTy>-<Hx,y>)/<Hx,y> = "
                       << (zz1-zz2)/zz2 << std::endl;

    BOOST_CHECK(zz1 != zero);
    BOOST_CHECK(zz2 != zero);
    BOOST_CHECK_CLOSE(zz1, zz2, tol);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testTangentLinear() {
  // Test  ||(hop(x+alpha*dx)-hop(x)) - hoptl(alpha*dx)|| < tol
  typedef ObsTestsFixture<MODEL>         Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::ObsAuxIncrement<MODEL>   ObsAuxIncr_;
  typedef oops::LinearObsOperator<MODEL> LinearObsOperator_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  const double tol = TestEnvironment::config().getDouble("LinearObsOpTest.toleranceTL");
  const int iter = TestEnvironment::config().getDouble("LinearObsOpTest.testiterTL");

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    ObsOperator_ hop(Test_::obspace()[jj], conf[jj]);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    LinearObsOperator_ hoptl(Test_::obspace()[jj], conf[jj]);

    // read geovals from the file
    const eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    const GeoVaLs_ x0(gconf, hop.variables());
    GeoVaLs_ x(gconf, hop.variables());

    // initialize obs bias
    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);
    const ObsAuxIncr_ ybinc(biasConf);

    // set TL trajectory to the geovals from the file
    hoptl.setTrajectory(x0, ybias);

    // create obsvectors
    ObsVector_ y1(Test_::obspace()[jj], hop.observed());
    ObsVector_ y2(Test_::obspace()[jj], hop.observed());
    ObsVector_ y3(Test_::obspace()[jj], hop.observed());

    // y1 = hop(x0)
    hop.simulateObs(x0, y1, ybias);

    // randomize dx
    GeoVaLs_ dx(gconf, hoptl.variables());
    dx.random();

    double alpha = 0.1;
    for (int jter = 0; jter < iter; ++jter) {
      // x = x0 + alpha*dx
      dx *= alpha;
      x = x0;
      x += dx;

      // y2 = hop(x0+alpha*dx)
      hop.simulateObs(x, y2, ybias);
      y2 -= y1;
      // y3 = hoptl(alpha*dx)
      hoptl.simulateObsTL(dx, y3, ybinc);
      y2 -= y3;
      double test_norm = y2.rms();
      oops::Log::info() << "Iter:" << jter << " ||(h(x+alpha*dx)-h(x))/h'(alpha*dx)||="
                         << test_norm << std::endl;
    }
    BOOST_CHECK(y2.rms() < tol);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> class LinearObsOperator : public oops::Test {
 public:
  LinearObsOperator() {}
  virtual ~LinearObsOperator() {}
 private:
  std::string testid() const {return "test::LinearObsOperator<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/LinearObsOperator");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testLinearity<MODEL>));
    ts->add(BOOST_TEST_CASE(&testTangentLinear<MODEL>));
    ts->add(BOOST_TEST_CASE(&testAdjoint<MODEL>));
    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LINEAROBSOPERATOR_H_
