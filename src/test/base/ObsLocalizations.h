/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_BASE_OBSLOCALIZATIONS_H_
#define TEST_BASE_OBSLOCALIZATIONS_H_

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Departures.h"
#include "oops/base/Geometry.h"
#include "oops/base/ObsLocalizations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

/// \brief Tests ObsLocalizations::computeLocalization method.
/// \details Tests that for obs localization around specified in yaml Geometry points:
/// 1. number of local obs matches reference ("reference local nobs")
/// 2. obs localization applied to a vector of ones computes a reference value for rmse
/// Reference gridpoints are specified in yaml as "observation localizations reference:"
/// includes: gridpoints.{lons, lats, depth}, rmse, local nobs
/// They don't have to be exactly equal to the lon/lat of the Geometry gridpoints,
/// but should be no further than 1.e-5 distance away.
/// test loops through all obs spaces and obs localizations
template <typename MODEL, typename OBS> void testObsLocalizations() {
  typedef ObsTestsFixture<OBS>                   Test_;
  typedef oops::Departures<OBS>                  Observations_;
  typedef oops::Geometry<MODEL>                  Geometry_;
  typedef oops::GeometryIterator<MODEL>          GeometryIterator_;
  typedef oops::ObsLocalizations<MODEL, OBS>     ObsLocalizations_;
  typedef oops::ObsSpaces<OBS>                   ObsSpaces_;

  // configure geometry
  const eckit::LocalConfiguration geometryConfig(TestEnvironment::config(), "geometry");
  Geometry_ geometry(geometryConfig, oops::mpi::world());

  // configure obs spaces and obs localizations
  const eckit::LocalConfiguration obsConfig(TestEnvironment::config(), "observations");
  const ObsSpaces_ & obspaces = Test_::obspace();
  ObsLocalizations_ obslocs(obsConfig, obspaces);
  oops::Log::test() << "Testing obs-space localizations: " << obslocs << std::endl;

  // get reference values
  const eckit::LocalConfiguration referenceConfig(TestEnvironment::config(),
                "observation localizations reference");
  const std::vector<double> lons = referenceConfig.getDoubleVector("gridpoints.lons");
  const std::vector<double> lats = referenceConfig.getDoubleVector("gridpoints.lats");
  const std::vector<double> depths = referenceConfig.getDoubleVector("gridpoints.depths");
  const std::vector<size_t> nobs_local_ref = referenceConfig.getUnsignedVector("local nobs");
  const std::vector<double> ref_locvector_rms = referenceConfig.getDoubleVector("rms");
  ASSERT(lons.size() == lats.size());
  ASSERT(lons.size() == depths.size());
  ASSERT(lons.size() == nobs_local_ref.size());
  ASSERT(lons.size() > 0);
  std::vector<eckit::geometry::Point3> reference_points;
  for (size_t jpoint = 0; jpoint < lons.size(); ++jpoint) {
    reference_points.emplace_back(lons[jpoint], lats[jpoint], depths[jpoint]);
  }

  Observations_ locvector(obspaces);
  Observations_ obsvector(obspaces);

  size_t total_tested = 0;
  std::vector<size_t> nobs_local(nobs_local_ref.size(), 0);
  std::vector<double> locvector_rms(nobs_local_ref.size(), 0);

  // loop over geometry points
  for (GeometryIterator_ ii = geometry.begin(); ii != geometry.end(); ++ii) {
    // debug print to help decide which points to specify for reference
    // set OOPS_DEBUG environment variable to -1 to see prints from all MPI tasks
    if (referenceConfig.getBool("print iterator", false)) {
      oops::Log::debug() << "Iterating over " << std::setprecision(9) << ii << ": "
                       << ii << std::endl;
    }
    // check if we need to test at this location (if there are any points in the
    // reference point list within 1e-5 of this locationn)
    const auto & it = std::find_if(reference_points.begin(), reference_points.end(),
          [ii] (const eckit::geometry::Point3 & point) {return point.distance(*ii) < 1e-5;});
    if ( it != reference_points.end() && (total_tested < lons.size()) ) {
      total_tested++;
      size_t index = it - reference_points.begin();
      locvector.ones();
      obslocs.computeLocalization(ii, locvector);
      oops::Log::test() << "Obs localization with geometry iterator: " << ii << ": "
                        << *ii << std::endl;
      oops::Log::test() << "Localization values: " << locvector << std::endl;
      oops::Log::test() << "Local vector nobs and reference: " << locvector.nobs() << ", "
                        << nobs_local_ref[index] << std::endl;
      oops::Log::debug() << "Local vector stats lat,lon,nobs,rms: " <<  *ii << ", "
                         << locvector.nobs() << ", " << locvector.rms() << std::endl;

      // save number of local obs to be tested later
      nobs_local[index] = locvector.nobs();

      // apply localization to a vector of ones
      obsvector.ones();
      obsvector *= locvector;
      oops::Log::test() << "Localization applied to local ObsVector of ones: " <<
                           obsvector << std::endl;
      // save localized vector rms to be tested later
      locvector_rms[index] = obsvector.rms();
    }
  }
  // gather number of tested observations, nobs_local and locvector_rms
  // (only one MPI task owns one gridpoint)
  Test_::comm().allReduceInPlace(total_tested, eckit::mpi::sum());
  Test_::comm().allReduceInPlace(nobs_local.begin(), nobs_local.end(), eckit::mpi::sum());
  Test_::comm().allReduceInPlace(locvector_rms.begin(), locvector_rms.end(), eckit::mpi::sum());
  // check that we tested all gridpoints
  EXPECT_EQUAL(total_tested, lons.size());
  // Test that computed number of local obs is the same as reference
  EXPECT_EQUAL(nobs_local_ref, nobs_local);
  // check value of the rms is close to reference
  ASSERT(lons.size() == ref_locvector_rms.size());
  oops::Log::debug() << "reference RMS" << ref_locvector_rms << std::endl;
  oops::Log::debug() << "computed RMS" << locvector_rms << std::endl;
  for (size_t jpoint = 0; jpoint < nobs_local.size(); ++jpoint) {
    if (nobs_local[jpoint] > 0) {
      EXPECT(std::abs(locvector_rms[jpoint]-ref_locvector_rms[jpoint]) < 1.e-5);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS> class ObsLocalizations : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;

 public:
  using oops::Test::Test;
  virtual ~ObsLocalizations() = default;
 private:
  std::string testid() const override {return "test::ObsLocalizations<" + MODEL::name() + ","
                                               + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("base/ObsLocalizations/testObsLocalizations")
      { testObsLocalizations<MODEL, OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_BASE_OBSLOCALIZATIONS_H_
