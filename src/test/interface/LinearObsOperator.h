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

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/scoped_ptr.hpp>

#include "eckit/testing/Test.h"
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
    eckit::LocalConfiguration obsopconf(conf[jj], "ObsOperator");
    boost::scoped_ptr<LinearObsOperator_> ov(
      new LinearObsOperator_(Test_::obspace()[jj], obsopconf));
    EXPECT(ov.get());

    ov.reset();
    EXPECT(!ov.get());
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
  const double tol = 1.0e-11;
  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    eckit::LocalConfiguration obsopconf(conf[jj], "ObsOperator");
    ObsOperator_ hop(Test_::obspace()[jj], obsopconf);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    LinearObsOperator_ hoptl(Test_::obspace()[jj], obsopconf);

    // read geovals from the file
    const eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    const GeoVaLs_ gval(gconf, hop.variables());

    // initialize obs bias
    const ObsAuxCtrl_ ybias(conf[jj]);
    const ObsAuxIncr_ ybinc(conf[jj]);

    // set trajectory for TL/AD to be the geovals from the file
    hoptl.setTrajectory(gval, ybias);

    // create obsvector
    ObsVector_ dy1(Test_::obspace()[jj], hop.observed());

    // create geovals
    GeoVaLs_ dx(gconf, hoptl.variables());

    // test rms(Hdx) = 0, when dx = 0
    dx.zero();
    hoptl.simulateObsTL(dx, dy1, ybinc);
    EXPECT(dy1.rms() == zero);

    // test rms(Hdx) > 0, when dx is random
    dx.random();
    hoptl.simulateObsTL(dx, dy1, ybinc);
    EXPECT(dy1.rms() > zero);

    // test k * H * dx ~ H * (k*dx)
    dy1 *= coef;
    dx  *= coef;
    ObsVector_ dy2(Test_::obspace()[jj], hop.observed());
    hoptl.simulateObsTL(dx, dy2, ybinc);

    dy1 -= dy2;

    EXPECT(dy1.rms()< tol);
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
    eckit::LocalConfiguration obsopconf(conf[jj], "ObsOperator");
    ObsOperator_ hop(Test_::obspace()[jj], obsopconf);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    LinearObsOperator_ hoptl(Test_::obspace()[jj], obsopconf);

    // read geovals from the file
    eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    const GeoVaLs_ gval(gconf, hop.variables());

    // initialize bias correction
    const ObsAuxCtrl_ ybias(conf[jj]);
    ObsAuxIncr_ ybinc(conf[jj]);

    // set TL/AD trajectory to the geovals from the file
    hoptl.setTrajectory(gval, ybias);

    ObsVector_ dy1(Test_::obspace()[jj], hop.observed());
    ObsVector_ dy2(Test_::obspace()[jj], hop.observed());
    GeoVaLs_ dx1(gconf, hoptl.variables());
    GeoVaLs_ dx2(gconf, hoptl.variables());

    // calculate dy1 = H dx1 (with random dx1)
    dx1.random();
    EXPECT(dot_product(dx1, dx1) > zero);  //  BOOST_REQUIRE
    hoptl.simulateObsTL(dx1, dy1, ybinc);
    EXPECT(dot_product(dy1, dy1) > zero);

    // calculate dx2 = HT dy2 (with random dy2)
    dy2.random();
    EXPECT(dot_product(dy2, dy2) > zero);  //  BOOST_REQUIRE
    dx2.zero();
    hoptl.simulateObsAD(dx2, dy2, ybinc);
    EXPECT(dot_product(dx2, dx2) > zero);

    const double zz1 = dot_product(dx1, dx2);
    const double zz2 = dot_product(dy1, dy2);

    oops::Log::info() << "Adjoint test result: (<x,HTy>-<Hx,y>)/<Hx,y> = "
                       << (zz1-zz2)/zz2 << std::endl;

    EXPECT(zz1 != zero);
    EXPECT(zz2 != zero);
    EXPECT(oops::is_close(zz1, zz2, tol));
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
  const double alpha = TestEnvironment::config().getDouble("LinearObsOpTest.coefTL", 0.1);
  const int iter = TestEnvironment::config().getInt("LinearObsOpTest.testiterTL", 1);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    // initialize observation operator (set variables requested from the model,
    // variables simulated by the observation operator, other init)
    eckit::LocalConfiguration obsopconf(conf[jj], "ObsOperator");
    ObsOperator_ hop(Test_::obspace()[jj], obsopconf);
    // initialize TL/AD observation operator (set model variables for Jacobian),
    // other init)
    LinearObsOperator_ hoptl(Test_::obspace()[jj], obsopconf);

    // read geovals from the file
    const eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    const GeoVaLs_ x0(gconf, hop.variables());
    GeoVaLs_ x(gconf, hop.variables());

    // initialize obs bias
    const ObsAuxCtrl_ ybias(conf[jj]);
    const ObsAuxIncr_ ybinc(conf[jj]);

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
      oops::Log::info() << "Iter:" << jter << " ||(h(x+alpha*dx)-h(x)-h'*(alpha*dx))||="
                        << test_norm << std::endl;
    }
    EXPECT(y2.rms() < tol);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearObsOperator : public oops::Test {
 public:
  LinearObsOperator() {}
  virtual ~LinearObsOperator() {}
 private:
  std::string testid() const {return "test::LinearObsOperator<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/GeometryIterator/testConstructor")
      { testConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/GeometryIterator/testLinearity")
      { testLinearity<MODEL>(); });
    ts.emplace_back(CASE("interface/GeometryIterator/testTangentLinear")
      { testTangentLinear<MODEL>(); });
    ts.emplace_back(CASE("interface/GeometryIterator/testAdjoint")
      { testAdjoint<MODEL>(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LINEAROBSOPERATOR_H_
