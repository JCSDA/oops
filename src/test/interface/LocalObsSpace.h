/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_LOCALOBSSPACE_H_
#define TEST_INTERFACE_LOCALOBSSPACE_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/geometry/Point2.h"
#include "eckit/testing/Test.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/runs/Test.h"
#include "oops/util/dot_product.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> void testLocalObsSpace() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::ObsSpace<MODEL>               LocalObsSpace_;
  typedef oops::ObsVector<MODEL>              ObsVector_;

  const eckit::LocalConfiguration localconf(TestEnvironment::config(), "LocalObsSpace");

  // get center (for localization) from yaml
  eckit::LocalConfiguration geolocconf(localconf, "GeoLocation");
  double lon = geolocconf.getDouble("lon");
  double lat = geolocconf.getDouble("lat");
  const eckit::geometry::Point2 center(lon, lat);

  // get localization info from yaml
  eckit::LocalConfiguration locconf(localconf, "Localization");

  int totalNobs = 0;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    // initialize local observation space
    LocalObsSpace_ localobs(Test_::obspace()[jj], center, locconf);
    oops::Log::info() << "Local obs within " << locconf << " from " << center <<
                         ": " << localobs << std::endl;

    // count local nobs
    ObsVector_ localvec(localobs);
    ObsVector_ globvec(Test_::obspace()[jj]);
    totalNobs += localvec.nobs();
  }

  // test that local nobs is equal to the reference value
  const int ref_nobs = localconf.getInt("reference nobs");
  EXPECT(totalNobs == ref_nobs);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testLocalObsVector() {
  typedef ObsTestsFixture<MODEL> Test_;
  typedef oops::ObsSpace<MODEL>               LocalObsSpace_;
  typedef oops::ObsVector<MODEL>              ObsVector_;

  const eckit::LocalConfiguration localconf(TestEnvironment::config(), "LocalObsSpace");

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    // get center (for localization) from yaml
    eckit::LocalConfiguration geolocconf(localconf, "GeoLocation");
    double lon = geolocconf.getDouble("lon");
    double lat = geolocconf.getDouble("lat");
    const eckit::geometry::Point2 center(lon, lat);

    // get localization info from yaml
    eckit::LocalConfiguration locconf(localconf, "Localization");

    const std::string varname = localconf.getString("varname");

    // initialize full ObsVector for a specified variable
    ObsVector_ fullvec(Test_::obspace()[jj], varname);
    oops::Log::info() << "Full Obsvector: " << fullvec << std::endl;

    // initialize local observation space
    LocalObsSpace_ localobs(Test_::obspace()[jj], center, locconf);
    oops::Log::info() << "Local obs within " << locconf << " from " << center <<
                         ": " << localobs << std::endl;

    // intialize local obsvector by reading specified variable from local obsspace
    ObsVector_ localvec1(localobs, varname);
    oops::Log::info() << "Local Obsvector from Local Obsspace: " << localvec1 << std::endl;
    // initialize local obsvector from full obsvector using local obsspace
    ObsVector_ localvec2(localobs, fullvec);
    oops::Log::info() << "Local ObsVector from full ObsVector: " << localvec2 << std::endl;
    // check that the two are equal
    EXPECT(localvec1.nobs() == localvec2.nobs());
    localvec2 -= localvec1;
    const double rms = dot_product(localvec2, localvec2);
    EXPECT(rms == 0);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL> class LocalObsSpace : public oops::Test {
 public:
  LocalObsSpace() {}
  virtual ~LocalObsSpace() {}
 private:
  std::string testid() const {return "test::LocalObsSpace<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();
    ts.emplace_back(CASE("interface/LocalObsSpace/testLocalObsSpace")
      { testLocalObsSpace<MODEL>(); });
    ts.emplace_back(CASE("interface/LocalObsSpace/testLocalObsVector")
      { testLocalObsVector<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_LOCALOBSSPACE_H_
