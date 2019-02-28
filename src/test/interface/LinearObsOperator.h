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

using eckit::types::is_approximately_equal;

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
    EXPECT(ov.get());

    ov.reset();
    EXPECT(!ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearity() {
  typedef ObsTestsFixture<MODEL>         Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::Locations<MODEL>         Locations_;
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
    ObsOperator_ hop(Test_::obspace()[jj], conf[jj]);
    LinearObsOperator_ hoptl(Test_::obspace()[jj], conf[jj]);

    const eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    Locations_ locs(hop.locations(Test_::tbgn(), Test_::tend()));
    const GeoVaLs_ gval(gconf, hop.variables());

    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);
    hoptl.setTrajectory(gval, ybias);

    const ObsAuxIncr_ ybinc(biasConf);
    ObsVector_ dy1(Test_::obspace()[jj], hop.observed());
    GeoVaLs_ gv(gconf, hoptl.variables());

    gv.zero();
    hoptl.simulateObsTL(gv, dy1, ybinc);

    EXPECT(dy1.rms() == zero);

    gv.random();
    hoptl.simulateObsTL(gv, dy1, ybinc);
    dy1 *= coef;
    EXPECT(dy1.rms() > zero);

    gv *= coef;
    ObsVector_ dy2(Test_::obspace()[jj], hop.observed());
    hoptl.simulateObsTL(gv, dy2, ybinc);

    dy1 -= dy2;

    EXPECT(dy1.rms() < tol);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testAdjoint() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::Locations<MODEL>         Locations_;
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
    ObsOperator_ hop(Test_::obspace()[jj], conf[jj]);
    LinearObsOperator_ hoptl(Test_::obspace()[jj], conf[jj]);
    eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    Locations_ locs(hop.locations(Test_::tbgn(), Test_::tend()));
    const GeoVaLs_ gval(gconf, hop.variables());

    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);

    hoptl.setTrajectory(gval, ybias);

    ObsAuxIncr_ ybinc(biasConf);

    ObsVector_ dy1(Test_::obspace()[jj], hop.observed());
    ObsVector_ dy2(Test_::obspace()[jj], hop.observed());
    GeoVaLs_ gv1(gconf, hoptl.variables());
    GeoVaLs_ gv2(gconf, hoptl.variables());

    gv1.random();
    EXPECT(dot_product(gv1, gv1) > zero);  // BOOST_REQUIRE
    hoptl.simulateObsTL(gv1, dy1, ybinc);
    EXPECT(dot_product(dy1, dy1) > zero);

    dy2.random();
    EXPECT(dot_product(dy2, dy2) > zero);  // BOOST_REQUIRE
    gv2.zero();
    hoptl.simulateObsAD(gv2, dy2, ybinc);
    EXPECT(dot_product(gv2, gv2) > zero);

    const double zz1 = dot_product(gv1, gv2);
    const double zz2 = dot_product(dy1, dy2);

    oops::Log::debug() << "Adjoint test result: (<x,HTy>-<Hx,y>)/<Hx,y> = "
                       << (zz1-zz2)/zz2 << std::endl;

    EXPECT(zz1 != zero);
    EXPECT(zz2 != zero);
    EXPECT(is_approximately_equal(zz1, zz2, tol));
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testTangentLinear() {
  // Test  ||(hop(x+alpha*dx)-hop(x)) - hoptl(alpha*dx)|| < tol
  typedef ObsTestsFixture<MODEL>         Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::Locations<MODEL>         Locations_;
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
    LinearObsOperator_ hoptl(Test_::obspace()[jj], conf[jj]);
    ObsOperator_ hop(Test_::obspace()[jj], conf[jj]);

    const eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    Locations_ locs(hop.locations(Test_::tbgn(), Test_::tend()));

    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);

    const ObsAuxIncr_ ybinc(biasConf);

    ObsVector_ y1(Test_::obspace()[jj], hop.observed());   // y1 = hop(x)
    ObsVector_ y2(Test_::obspace()[jj], hop.observed());   // y2 = hop(x+alpha*dx)
    ObsVector_ y3(Test_::obspace()[jj], hop.observed());   // y3 = hoptl(alpha*dx)

    GeoVaLs_ gv(gconf, hop.variables());  // Background

    hoptl.setTrajectory(gv, ybias);

    hop.simulateObs(gv, y1, ybias);

    GeoVaLs_ dgv(gconf, hoptl.variables());
    dgv.random();

    GeoVaLs_ gv0(gconf, hop.variables());
    gv0 = gv;
    ObsVector_ y3_init(Test_::obspace()[jj], hop.observed());
    y3_init = y3;
    double alpha = 0.1;
    for (int jter = 0; jter < iter; ++jter) {
      gv = gv0;
      dgv *= alpha;
      gv += dgv;

      hop.simulateObs(gv, y2, ybias);
      y2 -= y1;
      hoptl.simulateObsTL(dgv, y3, ybinc);
      y2 -= y3;
      double test_norm = y2.rms();
      y3 = y3_init;
      oops::Log::debug() << "Iter:" << jter << " ||(h(x+alpha*dx)-h(x))/h'(alpha*dx)||="
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
