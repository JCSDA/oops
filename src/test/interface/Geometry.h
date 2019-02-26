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

#include <memory>  // for std::unique_ptr
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>


#include "eckit/config/Configuration.h"
#include "eckit/testing/Test.h"
#include "oops/interface/Geometry.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"

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

  EXPECT(geom.get());

  geom.reset();
  EXPECT(!geom.get());
}

// -----------------------------------------------------------------------------
template <typename MODEL> void testCopyConstructor() {
  typedef oops::Geometry<MODEL>        Geometry_;
  boost::scoped_ptr<Geometry_> geom(new Geometry_(GeometryFixture<MODEL>::getConfig()));


  boost::scoped_ptr<Geometry_> other(new Geometry_(*geom));
  EXPECT(other.get());

  other.reset();
  EXPECT(!other.get());

  EXPECT(geom.get());
}

// -----------------------------------------------------------------------------
template <typename MODEL> class Geometry : public oops::Test {
 public:
  Geometry() {}
  virtual ~Geometry() {}
 private:
  std::string testid() const {return "test::Geometry<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/Geometry/testConstructor")
      { testConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/Geometry/testCopyConstructor")
      { testCopyConstructor<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GEOMETRY_H_

