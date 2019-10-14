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
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef oops::Locations<MODEL>        Locations_;

  const eckit::LocalConfiguration conf(TestEnvironment::config(), "Locations");
  std::unique_ptr<Locations_> locs(new Locations_(conf, oops::mpi::comm()));
  EXPECT(locs.get());

  locs.reset();
  EXPECT(!locs.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef oops::Locations<MODEL>        Locations_;

  const eckit::LocalConfiguration conf(TestEnvironment::config(), "Locations");
  std::unique_ptr<Locations_> locs(new Locations_(conf, oops::mpi::comm()));
  EXPECT(locs.get());

  std::unique_ptr<Locations_> other_locs(new Locations_(*locs));
  EXPECT(other_locs.get());

  locs.reset();
  EXPECT(!locs.get());

  other_locs.reset();
  EXPECT(!other_locs.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class Locations : public oops::Test {
 public:
  Locations() {}
  virtual ~Locations() {}
 private:
  std::string testid() const {return "test::Locations<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/Locations/testConstructor")
      { testConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/Locations/testCopyConstructor")
      { testCopyConstructor<MODEL>(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_LOCATIONS_H_
