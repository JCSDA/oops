/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_GEOMETRY_H_
#define TEST_INTERFACE_GEOMETRY_H_

#include <string>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/Geometry.h"
#include "test/TestEnvironment.h"
#include "eckit/config/Configuration.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class GeometryFixture : private boost::noncopyable {
 public:
  static const eckit::Configuration & getConfig() {return *getInstance().conf_;}

 private:
  static GeometryFixture<MODEL>& getInstance() {
    static GeometryFixture<MODEL> theGeometryFixture;
    return theGeometryFixture;
  }

  GeometryFixture() {
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "Geometry"));
  }

  ~GeometryFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
};
// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef oops::Geometry<MODEL>        Geometry_;

  boost::scoped_ptr<Geometry_> geom(new Geometry_(GeometryFixture<MODEL>::getConfig()));
  BOOST_CHECK(geom.get());

  geom.reset();
  BOOST_CHECK(!geom.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef oops::Geometry<MODEL>        Geometry_;
  boost::scoped_ptr<Geometry_> geom(new Geometry_(GeometryFixture<MODEL>::getConfig()));

  boost::scoped_ptr<Geometry_> other(new Geometry_(*geom));
  BOOST_CHECK(other.get());

  other.reset();
  BOOST_CHECK(!other.get());

  BOOST_CHECK(geom.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> class Geometry : public oops::Test {
 public:
  Geometry() {}
  virtual ~Geometry() {}
 private:
  std::string testid() const {return "test::Geometry<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/Geometry");

    ts->add(BOOST_TEST_CASE(&testConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testCopyConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GEOMETRY_H_
