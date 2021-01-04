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

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/Configuration.h"
#include "eckit/testing/Test.h"
#include "oops/interface/Geometry.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "test/interface/GeometryFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// \brief Tests constructor and print method
template <typename MODEL> void testConstructor() {
  typedef oops::Geometry<MODEL>        Geometry_;

  std::unique_ptr<Geometry_> geom(new Geometry_(GeometryFixture<MODEL>::getParameters(),
                                                oops::mpi::world(), oops::mpi::myself()));
  EXPECT(geom.get());
  oops::Log::test() << "Testing geometry: " << *geom << std::endl;
  geom.reset();
  EXPECT(!geom.get());
}

// -----------------------------------------------------------------------------
template <typename MODEL> class Geometry : public oops::Test {
 public:
  Geometry() {}
  virtual ~Geometry() {}
 private:
  std::string testid() const override {return "test::Geometry<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/Geometry/testConstructor")
      { testConstructor<MODEL>(); });
  }

  void clear() const override {}
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GEOMETRY_H_
