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

#include "oops/runs/Test.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/Variables.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"
#include "util/dot_product.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class GeoVaLsFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>          Geometry_;
  typedef oops::ObservationSpace<MODEL>  ObsSpace_;
  typedef oops::Variables<MODEL>         Variables_;

 public:
  static const ObsSpace_   & obspace() {return *getInstance().obspace_;}
  static const Variables_  & vars()    {return *getInstance().vars_;}
  static const Geometry_   & resol()   {return *getInstance().resol_;}
  static const util::DateTime & t1()   {return *getInstance().t1_;}
  static const util::DateTime & t2()   {return *getInstance().t2_;}

 private:
  static GeoVaLsFixture<MODEL>& getInstance() {
    static GeoVaLsFixture<MODEL> theGeoVaLsFixture;
    return theGeoVaLsFixture;
  }

  GeoVaLsFixture() {
    t1_.reset(new util::DateTime(TestEnvironment::config().getString("window_begin")));
    t2_.reset(new util::DateTime(TestEnvironment::config().getString("window_end")));

    const eckit::LocalConfiguration varConfig(TestEnvironment::config(), "Variables");
    vars_.reset(new Variables_(varConfig));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));

    std::vector<eckit::LocalConfiguration> obsConfs;
    TestEnvironment::config().get("Observations", obsConfs);
    BOOST_CHECK(obsConfs.size() > 0);
    const eckit::LocalConfiguration obsConf(obsConfs[0], "Observation");
    obspace_.reset(new ObsSpace_(obsConf, *t1_, *t2_));
  }

  ~GeoVaLsFixture() {}

  boost::scoped_ptr<ObsSpace_>        obspace_;
  boost::scoped_ptr<const Variables_> vars_;
  boost::scoped_ptr<const Geometry_>  resol_;
  boost::scoped_ptr<const util::DateTime> t1_;
  boost::scoped_ptr<const util::DateTime> t2_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef GeoVaLsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>  GeoVaLs_;

  boost::scoped_ptr<GeoVaLs_> ov(new GeoVaLs_(Test_::obspace(), Test_::vars(),
                                              Test_::t1(), Test_::t2(), Test_::resol()));
  BOOST_CHECK(ov.get());

  ov.reset();
  BOOST_CHECK(!ov.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testUtils() {
  typedef GeoVaLsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>  GeoVaLs_;

  GeoVaLs_ gval(Test_::obspace(), Test_::vars(), Test_::t1(), Test_::t2(), Test_::resol());

  gval.random();
  const double zz1 = dot_product(gval, gval);
  BOOST_CHECK(zz1 > 0.0);

  gval.zero();
  const double zz2 = dot_product(gval, gval);
  BOOST_CHECK_EQUAL(zz2, 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testRead() {
  typedef oops::GeoVaLs<MODEL>  GeoVaLs_;

  std::vector<eckit::LocalConfiguration> conf;
  TestEnvironment::config().get("GeoVaLs", conf);

  const double tol = 1.0e-8;
  for (std::size_t jj = 0; jj < conf.size(); ++jj) {
    GeoVaLs_ gval(conf[jj]);
    const double xx = conf[jj].getDouble("norm");
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
