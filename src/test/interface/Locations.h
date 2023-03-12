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

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/interface/Locations.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// \brief tests test constructor and print method
template <typename OBS> void testConstructor() {
  typedef oops::Locations<OBS>        Locations_;

  const eckit::LocalConfiguration conf(TestEnvironment::config(), "locations");
  std::unique_ptr<Locations_> locs(new Locations_(conf, oops::mpi::world()));
  EXPECT(locs.get());
  oops::Log::info() << "Testing locations: " << *locs << std::endl;
  locs.reset();
  EXPECT(!locs.get());
}

// -----------------------------------------------------------------------------

template <typename OBS>
class Locations : public oops::Test {
 public:
  Locations() {}
  virtual ~Locations() {}
 private:
  std::string testid() const override {return "test::Locations<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/Locations/testConstructor")
      { testConstructor<OBS>(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LOCATIONS_H_
