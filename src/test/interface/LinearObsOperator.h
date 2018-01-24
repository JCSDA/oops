/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_INTERFACE_LINEAROBSOPERATOR_H_
#define TEST_INTERFACE_LINEAROBSOPERATOR_H_

#include <string>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "test/TestEnvironment.h"
#include "test/interface/ObsTestsFixture.h"
#include "util/dot_product.h"

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

template <typename MODEL> void testAdjoint() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::LinearObsOperator<MODEL> LinearObsOperator_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::ObsAuxIncrement<MODEL>   ObsAuxIncr_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

  const double zero = 0.0;
  const double tol = 1.0e-12;
  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);
  
  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    LinearObsOperator_ hop(Test_::obspace()[jj]);
    oops::Variables vars = hop.variables();

    eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    gconf.set("Variables", vars.asConfig());

    const GeoVaLs_ gval(gconf);

    const eckit::LocalConfiguration biasconf(conf[jj], "ObsBias");
    const ObsAuxCtrl_ ybias(biasconf);

    hop.setTrajectory(gval, ybias);

    ObsAuxIncr_ ybinc(biasconf);

    ObsVector_ dy1(Test_::obspace()[jj]);
    ObsVector_ dy2(Test_::obspace()[jj]);
    GeoVaLs_ gv1(gconf);
    GeoVaLs_ gv2(gconf);

    gv1.random();
    BOOST_REQUIRE(dot_product(gv1, gv1) > zero);
    hop.obsEquivTL(gv1, dy1, ybinc);
    BOOST_CHECK(dot_product(dy1, dy1) > zero);

    dy2.random();
    BOOST_REQUIRE(dot_product(dy2, dy2) > zero);
    hop.obsEquivAD(gv2, dy2, ybinc);
    BOOST_CHECK(dot_product(gv2, gv2) > zero);

    const double zz1 = dot_product(gv1, gv2);
    const double zz2 = dot_product(dy1, dy2);
    BOOST_CHECK(zz1 != zero);
    BOOST_CHECK(zz2 != zero);
    BOOST_CHECK_CLOSE(zz1, zz2, tol);
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
//    ts->add(BOOST_TEST_CASE(&testLinearity<MODEL>));
//    ts->add(BOOST_TEST_CASE(&testTangentLinear<MODEL>));    
    ts->add(BOOST_TEST_CASE(&testAdjoint<MODEL>));
    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LINEAROBSOPERATOR_H_
