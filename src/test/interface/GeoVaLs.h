/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_INTERFACE_GEOVALS_H_
#define TEST_INTERFACE_GEOVALS_H_

#include <cmath>
#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/base/ObsSpaces.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/dot_product.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL>
class GeoVaLsFixture : private boost::noncopyable {
  typedef oops::ObsSpaces<MODEL>  ObsSpaces_;

 public:
  static const util::DateTime  & tbgn() {return *getInstance().tbgn_;}
  static const util::DateTime  & tend() {return *getInstance().tend_;}
  static ObsSpaces_         & obspace() {return *getInstance().ospaces_;}
  static const oops::Variables & vars() {return *getInstance().vars_;}

 private:
  static GeoVaLsFixture<MODEL>& getInstance() {
    static GeoVaLsFixture<MODEL> theGeoVaLsFixture;
    return theGeoVaLsFixture;
  }

  GeoVaLsFixture(): tbgn_(), tend_(), ospaces_() {
    tbgn_.reset(new util::DateTime(TestEnvironment::config().getString("window_begin")));
    tend_.reset(new util::DateTime(TestEnvironment::config().getString("window_end")));

    const eckit::LocalConfiguration conf(TestEnvironment::config(), "Observations");
    ospaces_.reset(new ObsSpaces_(conf, *tbgn_, *tend_));

    const eckit::LocalConfiguration varConfig(TestEnvironment::config(), "Variables");
    vars_.reset(new oops::Variables(varConfig));
  }

  ~GeoVaLsFixture() {}

  boost::scoped_ptr<const util::DateTime> tbgn_;
  boost::scoped_ptr<const util::DateTime> tend_;
  boost::scoped_ptr<ObsSpaces_> ospaces_;
  boost::scoped_ptr<const oops::Variables> vars_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef GeoVaLsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>  GeoVaLs_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    boost::scoped_ptr<GeoVaLs_> ov(new GeoVaLs_(Test_::obspace()[jj], Test_::vars(),
                                                Test_::tbgn(), Test_::tend()));
    BOOST_CHECK(ov.get());

    ov.reset();
    BOOST_CHECK(!ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testUtils() {
  typedef GeoVaLsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>  GeoVaLs_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    GeoVaLs_ gval(Test_::obspace()[jj], Test_::vars(),
                  Test_::tbgn(), Test_::tend());

    gval.random();
    const double zz1 = dot_product(gval, gval);
    BOOST_CHECK(zz1 > 0.0);

    gval.zero();
    const double zz2 = dot_product(gval, gval);
    BOOST_CHECK_EQUAL(zz2, 0.0);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testRead() {
  typedef oops::GeoVaLs<MODEL>  GeoVaLs_;

  const eckit::LocalConfiguration obsconf(TestEnvironment::config(), "Observations");
  std::vector<eckit::LocalConfiguration> conf;
  obsconf.get("ObsTypes", conf);

  const double tol = 1.0e-8;
  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    const eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
    GeoVaLs_ gval(gconf);
    const double xx = gconf.getDouble("norm");
    const double zz = sqrt(dot_product(gval, gval));
    BOOST_CHECK_CLOSE(xx, zz, tol);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> class GeoVaLs : public oops::Test {
 public:
  GeoVaLs() {}
  virtual ~GeoVaLs() {}
 private:
  std::string testid() const {return "test::GeoVaLs<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/GeoVaLs");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testUtils<MODEL>));
    ts->add(BOOST_TEST_CASE(&testRead<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GEOVALS_H_
