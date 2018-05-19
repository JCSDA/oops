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

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    boost::scoped_ptr<LinearObsOperator_> ov(new LinearObsOperator_(Test_::obspace()[jj]));
    BOOST_CHECK(ov.get());

    ov.reset();
    BOOST_CHECK(!ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLinearity() {
  typedef ObsTestsFixture<MODEL>         Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::Locations<MODEL>         Locations_;
  typedef oops::LinearObsOperator<MODEL> LinearObsOperator_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::ObsAuxIncrement<MODEL>   ObsAuxIncr_;
  typedef oops::LinearObsOperator<MODEL> LinearObsOperator_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

  const double zero = 0.0;
  const double coef = 3.14;
  const double tol = TestEnvironment::config().getDouble("LinearObsOpTest.toleranceAD");
  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    LinearObsOperator_ hoptl(Test_::obspace()[jj]);

    const eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    Locations_ locs(Test_::obspace()[jj].locations(Test_::tbgn(), Test_::tend()));
    const GeoVaLs_ gval(locs, hoptl.variables(), gconf);

    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);
    hoptl.setTrajectory(gval, ybias);

    const ObsAuxIncr_ ybinc(biasConf);
    ObsVector_ dy1(Test_::obspace()[jj]);
    GeoVaLs_ gv(locs, hoptl.variables(), gconf);

    gv.zero();
    hoptl.obsEquivTL(gv, dy1, ybinc);

    BOOST_CHECK_EQUAL(dy1.rms(), zero);

    gv.random();
    hoptl.obsEquivTL(gv, dy1, ybinc);
    dy1 *= coef;
    BOOST_CHECK(dy1.rms() > zero);

    gv *= coef;
    ObsVector_ dy2(Test_::obspace()[jj]);
    hoptl.obsEquivTL(gv, dy2, ybinc);

    dy1 -= dy2;

    BOOST_CHECK_SMALL(dy1.rms(), tol);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testAdjoint() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::Locations<MODEL>         Locations_;
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
    LinearObsOperator_ hoptl(Test_::obspace()[jj]);
    eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    Locations_ locs(Test_::obspace()[jj].locations(Test_::tbgn(), Test_::tend()));
    const GeoVaLs_ gval(locs, hoptl.variables(), gconf);

    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);

    hoptl.setTrajectory(gval, ybias);

    ObsAuxIncr_ ybinc(biasConf);

    ObsVector_ dy1(Test_::obspace()[jj]);
    ObsVector_ dy2(Test_::obspace()[jj]);
    GeoVaLs_ gv1(locs, hoptl.variables(), gconf);
    GeoVaLs_ gv2(locs, hoptl.variables(), gconf);

    gv1.random();
    BOOST_REQUIRE(dot_product(gv1, gv1) > zero);
    hoptl.obsEquivTL(gv1, dy1, ybinc);
    BOOST_CHECK(dot_product(dy1, dy1) > zero);

    dy2.random();
    BOOST_REQUIRE(dot_product(dy2, dy2) > zero);
    hoptl.obsEquivAD(gv2, dy2, ybinc);
    BOOST_CHECK(dot_product(gv2, gv2) > zero);

    const double zz1 = dot_product(gv1, gv2);
    const double zz2 = dot_product(dy1, dy2);
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
  typedef oops::Locations<MODEL>         Locations_;
  typedef oops::LinearObsOperator<MODEL> LinearObsOperator_;
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
    LinearObsOperator_ hoptl(Test_::obspace()[jj]);
    ObsOperator_ hop(Test_::obspace()[jj]);

    const eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    Locations_ locs(Test_::obspace()[jj].locations(Test_::tbgn(), Test_::tend()));

    eckit::LocalConfiguration biasConf;
    conf[jj].get("ObsBias", biasConf);
    const ObsAuxCtrl_ ybias(biasConf);

    const ObsAuxIncr_ ybinc(biasConf);

    ObsVector_ y1(Test_::obspace()[jj]);   // y1 = hop(x)
    ObsVector_ y2(Test_::obspace()[jj]);   // y2 = hop(x+alpha*dx)
    ObsVector_ y3(Test_::obspace()[jj]);   // y3 = hoptl(alpha*dx)

    GeoVaLs_ gv(locs, hop.variables(), gconf);  // Background

    hoptl.setTrajectory(gv, ybias);

    hop.obsEquiv(gv, y1, ybias);

    GeoVaLs_ dgv(locs, hoptl.variables(), gconf);
    dgv.random();

    GeoVaLs_ gv0(locs, hop.variables(), gconf);
    gv0 = gv;
    ObsVector_ y3_init(Test_::obspace()[jj]);
    y3_init = y3;
    double alpha = 0.1;
    for (int jter = 0; jter < iter; ++jter) {
      gv = gv0;
      dgv *= alpha;
      gv += dgv;

      hop.obsEquiv(gv, y2, ybias);
      y2 -= y1;
      hoptl.obsEquivTL(dgv, y3, ybinc);
      y2 -= y3;
      double test_norm = y2.rms();
      y3 = y3_init;
      oops::Log::debug() << "Iter:" << jter << " ||(h(x+alpha*dx)-h(x))/h'(alpha*dx)||="
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
