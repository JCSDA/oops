/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_GEOMETRYITERATOR_H_
#define TEST_INTERFACE_GEOMETRYITERATOR_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/Configuration.h"
#include "eckit/testing/Test.h"

#include "oops/base/LocalIncrement.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/Increment.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"

#include "test/interface/Increment.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef oops::GeometryIterator<MODEL>   GeometryIterator_;
  typedef oops::Geometry<MODEL>           Geometry_;

  const eckit::LocalConfiguration geomConfig(TestEnvironment::config(), "geometry");
  Geometry_ geom(geomConfig, oops::mpi::world());

  std::unique_ptr<GeometryIterator_> iter(new GeometryIterator_(geom.begin()));
  EXPECT(iter.get());

  iter.reset();
  EXPECT(!iter.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testGetSetLocal() {
  typedef oops::Geometry<MODEL>          Geometry_;
  typedef oops::GeometryIterator<MODEL>  GeometryIterator_;
  typedef oops::Increment<MODEL>         Increment_;
  typedef IncrementFixture<MODEL>        Test_;

  const eckit::LocalConfiguration geomConfig(TestEnvironment::config(), "geometry");
  Geometry_ geom(geomConfig, oops::mpi::world());

  // randomize increment dx1
  Increment_ dx1(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx1.random();
  oops::Log::info() << "Increment dx1 (random): " << dx1 << std::endl;
  // zero out increment dx2
  Increment_ dx2(Test_::resol(), Test_::ctlvars(), Test_::time());
  dx2.zero();
  oops::Log::info() << "Increment dx2 (zero): " << dx2 << std::endl;

  for (GeometryIterator_ i = geom.begin(); i != geom.end(); ++i) {
    // get value for i-th gridpoint from dx1
    oops::LocalIncrement gp = dx1.getLocal(i);
    oops::Log::debug() << *i << gp << std::endl;
    // set this value for i-th gridpoint in dx2
    dx2.setLocal(gp, i);
  }
  oops::Log::info() << "Increment dx2 after dx2=dx1 (at every point): " << dx2 << std::endl;
  // compare two increments
  dx2 -= dx1;
  EXPECT(dx2.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> class GeometryIterator : public oops::Test {
 public:
  GeometryIterator() {}
  virtual ~GeometryIterator() {}

 private:
  std::string testid() const {return "test::GeometryIterator<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/GeometryIterator/testConstructor")
      { testConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/GeometryIterator/testGetSetLocal")
      { testGetSetLocal<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GEOMETRYITERATOR_H_
