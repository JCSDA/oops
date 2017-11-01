/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_LOCATIONS_H_
#define TEST_INTERFACE_LOCATIONS_H_

#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ObservationSpace.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class LocationsFixture : private boost::noncopyable {
  typedef oops::ObservationSpace<MODEL> ObsSpace_;
 public:
  static const ObsSpace_ & getObSpace() {return *getInstance().os_;}
  static const util::DateTime & gett1() {return *getInstance().t1_;}
  static const util::DateTime & gett2() {return *getInstance().t2_;}

 private:
  static LocationsFixture<MODEL>& getInstance() {
    static LocationsFixture<MODEL> theLocationsFixture;
    return theLocationsFixture;
  }

  LocationsFixture() {
    const util::DateTime bgn(TestEnvironment::config().getString("window_begin"));
    const util::DateTime end(TestEnvironment::config().getString("window_end"));

    std::vector<eckit::LocalConfiguration> obsConfs;
    TestEnvironment::config().get("Observations", obsConfs);
    BOOST_CHECK(obsConfs.size() > 0);
    const eckit::LocalConfiguration obsConf(obsConfs[0], "Observation");
    os_.reset(new ObsSpace_(obsConf, bgn, end));

    t1_.reset(new util::DateTime("2010-01-01T12:00:00Z"));
    t2_.reset(new util::DateTime("2010-01-02T00:00:00Z"));
  }

  ~LocationsFixture() {}

  boost::scoped_ptr<ObsSpace_> os_;
  boost::scoped_ptr<util::DateTime> t1_;
  boost::scoped_ptr<util::DateTime> t2_;
};
// -----------------------------------------------------------------------------

// template <typename MODEL> void testLocationsConstructor<MODEL> () {
// BOOST_TEST_CASE_TEMPLATE_FUNCTION( testLocationsConstructor, MODEL ) {
// BOOST_AUTO_TEST_CASE_TEMPLATE( testLocationsConstructor, MODEL, test_types ) {
template <typename MODEL> void testConstructors() {
  typedef oops::Locations<MODEL>        Locations_;

  boost::scoped_ptr<Locations_> locs(new Locations_(LocationsFixture<MODEL>::getObSpace(),
                                                    LocationsFixture<MODEL>::gett1(),
                                                    LocationsFixture<MODEL>::gett2()));

  BOOST_CHECK(locs.get());

  locs.reset();

  BOOST_CHECK(!locs.get());
}

// =============================================================================

template <typename MODEL> class Locations : public oops::Test {
 public:
  Locations() {}
  virtual ~Locations() {}
 private:
  std::string testid() const {return "test::Locations<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/Locations");

    ts->add(BOOST_TEST_CASE(&testConstructors<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_LOCATIONS_H_
