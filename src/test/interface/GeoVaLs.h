/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_INTERFACE_GEOVALS_H_
#define TEST_INTERFACE_GEOVALS_H_

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

 private:
  static GeoVaLsFixture<MODEL>& getInstance() {
    static GeoVaLsFixture<MODEL> theGeoVaLsFixture;
    return theGeoVaLsFixture;
  }

  GeoVaLsFixture() {
    const util::DateTime tbgn(TestEnvironment::config().getString("window_begin"));
    const util::DateTime tend(TestEnvironment::config().getString("window_end"));

    const eckit::LocalConfiguration varConfig(TestEnvironment::config(), "Variables");
    vars_.reset(new Variables_(varConfig));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));

    std::vector<eckit::LocalConfiguration> obsConfs;
    TestEnvironment::config().get("Observations", obsConfs);
    BOOST_CHECK(obsConfs.size() > 0);
    const eckit::LocalConfiguration obsConf(obsConfs[0], "Observation");
    obspace_.reset(new ObsSpace_(obsConf, tbgn, tend));
  }

  ~GeoVaLsFixture() {}

  boost::scoped_ptr<ObsSpace_>        obspace_;
  boost::scoped_ptr<const Variables_> vars_;
  boost::scoped_ptr<const Geometry_>  resol_;
};
// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef GeoVaLsFixture<MODEL> Test_;
  typedef oops::GeoVaLs<MODEL>  GeoVaLs_;

  const util::DateTime t1(TestEnvironment::config().getString("window_begin"));
  const util::DateTime t2(TestEnvironment::config().getString("window_end"));

  boost::scoped_ptr<GeoVaLs_> ov(new GeoVaLs_(Test_::obspace(), Test_::vars(),
                                              t1, t2, Test_::resol()));
  BOOST_CHECK(ov.get());

  ov.reset();
  BOOST_CHECK(!ov.get());
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

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GEOVALS_H_
