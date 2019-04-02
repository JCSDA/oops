/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TEST_INTERFACE_LOCALOBSERVATIONSPACE_H_
#define TEST_INTERFACE_LOCALOBSERVATIONSPACE_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/GeoDistance.h"
#include "oops/base/GeoLocation.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsVector.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testLocal() {
  typedef ObsVecFixture<MODEL> Test_;
  typedef oops::ObservationSpace<MODEL>       LocalObsSpace_;
  typedef oops::ObsVector<MODEL>              ObsVector_;

  const eckit::LocalConfiguration localconf(TestEnvironment::config(), "LocalObservationSpace");

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration geolocconf(localconf, "GeoLocation");
    const oops::GeoLocation center(geolocconf);
    eckit::LocalConfiguration distconf(localconf, "GeoDistance");
    const oops::GeoDistance dist(distconf);

    ObsVector_ fullvec(*Test_::obspace()[jj], Test_::observed(jj));
    oops::Log::info() << "Full Obsvector: " << fullvec << std::endl;

    LocalObsSpace_ local(*Test_::obspace()[jj], center, dist, -1);
    oops::Log::info() << "Local obs within " << dist << " from " << center <<
                         ": " <<local << std::endl;

    ObsVector_ localvec(local, Test_::observed(jj));
    oops::Log::info() << "Local Obsvector: " << localvec << std::endl;
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> class LocalObservationSpace : public oops::Test {
 public:
  LocalObservationSpace() {}
  virtual ~LocalObservationSpace() {}
 private:
  std::string testid() const {return "test::LocalObservationSpace<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/LocalObservationSpace/testLocal")
      { testLocal<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_LOCALOBSERVATIONSPACE_H_
