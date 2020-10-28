/*
 * (C) Copyright 2017-2020 UCAR
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
/// Tests that number of local observations in LocalObsSpace is the same as in reference
template <typename OBS> void testLocalObsSpace() {
  typedef ObsTestsFixture<OBS>         Test_;
  typedef oops::ObsSpace<OBS>          LocalObsSpace_;
  typedef oops::ObsVector<OBS>         ObsVector_;

  const eckit::LocalConfiguration localconf(TestEnvironment::config(), "local obs space");

  // get center (for localization) from yaml
  eckit::LocalConfiguration geolocconf(localconf, "location");
  double lon = geolocconf.getDouble("lon");
  double lat = geolocconf.getDouble("lat");
  const eckit::geometry::Point2 center(lon, lat);

  // get localization test from yaml
  eckit::LocalConfiguration locconf(localconf, "localization");

  // count local nobs for all obs types
  int totalNobs = 0;
  for (size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    // initialize local observation space
    LocalObsSpace_ localobs(Test_::obspace()[jj], center, locconf);
    oops::Log::test() << "Local obs within " << locconf << " from " << center <<
                         ": " << localobs << std::endl;

    // count local nobs
    ObsVector_ localvec(localobs);
    totalNobs += localvec.nobs();
  }

  // test that local nobs is equal to the reference value
  const int ref_nobs = localconf.getInt("reference nobs");
  EXPECT(totalNobs == ref_nobs);
}

// -----------------------------------------------------------------------------
/// Tests that constructing local ObsVector from local ObsSpace by reading
/// (ObsVector(const ObsSpace &, const std::string &) ctor), and by subsetting
/// full ObsVector (ObsVector(const ObsSpace &, const ObsVector &) ctor)
/// give the same results
template <typename OBS> void testLocalObsVector() {
  typedef ObsTestsFixture<OBS>              Test_;
  typedef oops::ObsSpace<OBS>               LocalObsSpace_;
  typedef oops::ObsVector<OBS>              ObsVector_;

  const eckit::LocalConfiguration localconf(TestEnvironment::config(), "local obs space");
  // get center (for localization) from yaml
  eckit::LocalConfiguration geolocconf(localconf, "location");
  double lon = geolocconf.getDouble("lon");
  double lat = geolocconf.getDouble("lat");
  const eckit::geometry::Point2 center(lon, lat);

  // get localization test from yaml
  eckit::LocalConfiguration locconf(localconf, "localization");
  const std::string varname = localconf.getString("variable name");

  for (size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    // initialize full ObsVector for a specified variable
    ObsVector_ fullvec(Test_::obspace()[jj], varname);
    oops::Log::test() << "Full Obsvector: " << fullvec << std::endl;

    // initialize local observation space
    LocalObsSpace_ localobs(Test_::obspace()[jj], center, locconf);
    oops::Log::test() << "Local obs within " << locconf << " from " << center <<
                           ": " << localobs << std::endl;

    // intialize local obsvector by reading specified variable from local obsspace
    ObsVector_ localvec1(localobs, varname);
    oops::Log::test() << "Local Obsvector from Local Obsspace: " << localvec1 << std::endl;

    // initialize local obsvector from full obsvector using local obsspace
    ObsVector_ localvec2(localobs, fullvec);
    oops::Log::test() << "Local ObsVector from full ObsVector: " << localvec2 << std::endl;
    // check that the two are equal
    EXPECT(localvec1.nobs() == localvec2.nobs());
    localvec2 -= localvec1;
    const double rms = dot_product(localvec2, localvec2);
    EXPECT(rms == 0);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> class LocalObsSpace : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;
 public:
  LocalObsSpace() {}
  virtual ~LocalObsSpace() {}
 private:
  std::string testid() const override {return "test::LocalObsSpace<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();
    ts.emplace_back(CASE("interface/LocalObsSpace/testLocalObsSpace")
      { testLocalObsSpace<OBS>(); });
    ts.emplace_back(CASE("interface/LocalObsSpace/testLocalObsVector")
      { testLocalObsVector<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_LOCALOBSSPACE_H_
