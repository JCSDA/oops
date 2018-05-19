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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/interface/Locations.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef oops::Locations<MODEL>        Locations_;

  const eckit::LocalConfiguration conf(TestEnvironment::config(), "Locations");
  boost::scoped_ptr<Locations_> locs(new Locations_(conf));
  BOOST_CHECK(locs.get());

  locs.reset();
  BOOST_CHECK(!locs.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef oops::Locations<MODEL>        Locations_;

  const eckit::LocalConfiguration conf(TestEnvironment::config(), "Locations");
  boost::scoped_ptr<Locations_> locs(new Locations_(conf));
  BOOST_CHECK(locs.get());

  boost::scoped_ptr<Locations_> other_locs(new Locations_(*locs));
  BOOST_CHECK(other_locs.get());

  locs.reset();
  BOOST_CHECK(!locs.get());

  other_locs.reset();
  BOOST_CHECK(!other_locs.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> class Locations : public oops::Test {
 public:
  Locations() {}
  virtual ~Locations() {}
 private:
  std::string testid() const {return "test::Locations<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/Locations");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testCopyConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LOCATIONS_H_
