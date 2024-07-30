/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_OBSLOCALIZATION_H_
#define TEST_INTERFACE_OBSLOCALIZATION_H_

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Geometry.h"
#include "oops/base/ObsLocalizationBase.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

/// \brief Tests ObsLocalization::computeLocalization method.
/// \details Tests that for obs localization around specified in yaml Geometry points:
/// 1. number of local obs matches reference ("reference local nobs")
/// 2. obs localization applied to a vector of ones computes a reference value for rmse
/// Reference gridpoints are specified in yaml as "reference gridpoints.{lons, lats, depth}"
/// They don't have to be exactly equal to the lon/lat of the Geometry gridpoints,
/// but should be no further than 1.e-5 distance away.
/// test loops through all obs spaces and obs localizations
template <typename MODEL, typename OBS> void testObsLocalization() {
  typedef ObsTestsFixture<OBS>                   Test_;
  typedef oops::Geometry<MODEL>                  Geometry_;
  typedef oops::GeometryIterator<MODEL>          GeometryIterator_;
  typedef oops::ObsLocalizationBase<MODEL, OBS>  ObsLocalization_;
  typedef oops::ObsSpace<OBS>                    ObsSpace_;
  typedef oops::ObsVector<OBS>                   ObsVector_;

  const eckit::LocalConfiguration geometryConfig(TestEnvironment::config(), "geometry");
  Geometry_ geometry(geometryConfig, oops::mpi::world());

  // loop over all obs spaces
  for (size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const ObsSpace_ & obspace = Test_::obspace()[jj];

    // retrieve and loop through all localizations in this obspace
    std::vector<eckit::LocalConfiguration> obsLocConfigs =
                        Test_::config(jj).getSubConfigurations("obs localizations");

    for (size_t oli = 0; oli < obsLocConfigs.size(); ++oli) {
      // initialize obs-space localization
      eckit::LocalConfiguration locconf(obsLocConfigs[oli]);

      // read reference local nobs values and reference gridpoints
      const std::vector<double> lons = locconf.getDoubleVector("reference gridpoints.lons");
      const std::vector<double> lats = locconf.getDoubleVector("reference gridpoints.lats");
      const std::vector<double> depths = locconf.getDoubleVector("reference gridpoints.depths");
      const std::vector<size_t> nobs_local_ref = locconf.getUnsignedVector("reference local nobs");
      ASSERT(lons.size() == lats.size());
      ASSERT(lons.size() == depths.size());
      ASSERT(lons.size() == nobs_local_ref.size());
      ASSERT(lons.size() > 0);
      std::vector<eckit::geometry::Point3> reference_points;
      for (size_t jpoint = 0; jpoint < lons.size(); ++jpoint) {
        reference_points.emplace_back(lons[jpoint], lats[jpoint], depths[jpoint]);
      }
      eckit::LocalConfiguration testconf;
      testconf.set("localization method", locconf.getString("localization method"));
      if (locconf.has("lengthscale"))
        testconf.set("lengthscale", locconf.getDouble("lengthscale"));
      if (locconf.has("soar horizontal decay"))
        testconf.set("soar horizontal decay", locconf.getDouble("soar horizontal decay"));

      if (locconf.has("apply log transformation"))
        testconf.set("apply log transformation", locconf.getBool("apply log transformation"));
      if (locconf.has("vertical lengthscale"))
        testconf.set("vertical lengthscale", locconf.getInt("vertical lengthscale"));
      if (locconf.has("ioda vertical coordinate"))
        testconf.set("ioda vertical coordinate", locconf.getString("ioda vertical coordinate"));
      if (locconf.has("ioda vertical coordinate group"))
        testconf.set("ioda vertical coordinate group",
                     locconf.getString("ioda vertical coordinate group"));
      if (locconf.has("localization function"))
        testconf.set("localization function", locconf.getString("localization function"));
      if (locconf.has("assign constant vertical coordinate to obs"))
        testconf.set("assign constant vertical coordinate to obs",
                     locconf.getBool("assign constant vertical coordinate to obs"));
      if (locconf.has("constant vertical coordinate value"))
        testconf.set("constant vertical coordinate value",
                     locconf.getDouble("constant vertical coordinate value"));

      if (locconf.has("base value"))
        testconf.set("base value", locconf.getDouble("base value"));
      if (locconf.has("rossby mult"))
        testconf.set("rossby mult", locconf.getDouble("rossby mult"));
      if (locconf.has("min grid mult"))
        testconf.set("min grid mult", locconf.getDouble("min grid mult"));
      if (locconf.has("min value"))
        testconf.set("min value", locconf.getDouble("min value"));
      if (locconf.has("max value"))
        testconf.set("max value", locconf.getDouble("max value"));
      if (locconf.has("search method"))
        testconf.set("search method", locconf.getString("search method"));
      if (locconf.has("distance type"))
        testconf.set("distance type", locconf.getString("distance type"));

      std::unique_ptr<ObsLocalization_> obsloc =
        oops::ObsLocalizationFactory<MODEL, OBS>::create(testconf, obspace);
      oops::Log::info() << "Testing obs-space localization: " << *obsloc << std::endl;

      ObsVector_ locvector(obspace);
      ObsVector_ obsvector(obspace);

      size_t total_tested = 0;
      std::vector<size_t> nobs_local(nobs_local_ref.size(), 0);
      std::vector<double> locvector_rms(nobs_local_ref.size(), 0);
      // loop over geometry points
      for (GeometryIterator_ ii = geometry.begin(); ii != geometry.end(); ++ii) {
        // check if we need to test at this location (if there are any points in the
        // reference point list within 1e-5 of this locationn)
        const auto & it = std::find_if(reference_points.begin(), reference_points.end(),
              [ii] (const eckit::geometry::Point3 & point) {return point.distance(*ii) < 1e-5;});
        // if iterator it within the bounds of the reference array, go ahead with further testing
        if ( it != reference_points.end() && (total_tested < lons.size()) ) {
          total_tested++;
          size_t index = it - reference_points.begin();
          locvector.ones();
          obsloc->computeLocalization(ii, locvector);
          oops::Log::info() << "Localization values: " << locvector << std::endl;
          oops::Log::debug() << "Local vector stats lat,lon,nobs,rms: " <<  *ii << ", "
                             << locvector.nobs() << ", " << locvector.rms() << std::endl;

          // save number of local obs to be tested later
          nobs_local[index] = locvector.nobs();

          // apply localization to a vector of ones
          obsvector.ones();
          obsvector *= locvector;
          oops::Log::info() << "Localization applied to local ObsVector of ones: " <<
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
      const std::vector<double> ref_locvector_rms = locconf.getDoubleVector("reference rms");
      ASSERT(lons.size() == ref_locvector_rms.size());
      oops::Log::debug() << "reference RMS" << ref_locvector_rms << std::endl;
      oops::Log::debug() << "computed RMS" << locvector_rms << std::endl;
      for (size_t jpoint = 0; jpoint < nobs_local.size(); ++jpoint) {
        if (nobs_local[jpoint] > 0) {
          EXPECT(std::abs(locvector_rms[jpoint]-ref_locvector_rms[jpoint]) < 1.e-5);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS> class ObsLocalization : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;

 public:
  using oops::Test::Test;
  virtual ~ObsLocalization() = default;

 private:
  std::string testid() const override {return "test::ObsLocalization<" + MODEL::name() + ","
                                               + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsLocalization/testObsLocalization")
      { testObsLocalization<MODEL, OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_OBSLOCALIZATION_H_
