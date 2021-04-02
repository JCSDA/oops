/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_OBSLOCALIZATION_H_
#define TEST_INTERFACE_OBSLOCALIZATION_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsLocalizationBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

/// Tests that obs localization applied to a zero vector returns zero.
/// Tests that obs localization applied to a vector of ones makes rms(obsvec) < 1
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
    eckit::LocalConfiguration locconf(Test_::config(jj), "obs localization");
    // variable used to check if localization is tested at least once below
    bool tested = false;
    // loop over geometry points
    for (GeometryIterator_ ii = geometry.begin(); ii != geometry.end(); ++ii) {
      // initialize local observation space and a local obs vector
      ObsSpace_ localobs(Test_::obspace()[jj], *ii, locconf);
      ObsVector_ obsvector(localobs);
      // only test if there are some observations in the local obs space
      if (obsvector.nobs() > 0) {
        tested = true;
        // initialize obs-space localization
        std::unique_ptr<ObsLocalization_> obsloc =
                        oops::ObsLocalizationFactory<MODEL, OBS>::create(locconf, localobs);
        oops::Log::test() << "Testing obs-space localization: " << *obsloc <<
                             " at geometry iterator " << ii << std::endl;
        // apply obs localization to a zero vector, check that result is zero
        obsvector.zero();
        EXPECT_EQUAL(obsvector.rms(), 0.0);
        ObsVector_ locvector(localobs);
        obsloc->computeLocalization(ii, locvector);
        obsvector *= locvector;
        EXPECT_EQUAL(obsvector.rms(), 0.0);
        // apply localization to a vector of ones, check that rms(result) < 1
        obsvector.ones();
        EXPECT_EQUAL(obsvector.rms(), 1.0);
        locvector.ones();
        obsloc->computeLocalization(ii, locvector);
        obsvector *= locvector;
        oops::Log::test() << "Localization applied to local ObsVector of ones: " <<
                             obsvector << std::endl;
        EXPECT(obsvector.rms() < 1.0);
      }
    }
    EXPECT(tested);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS> class ObsLocalization : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;

 public:
  ObsLocalization() = default;
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
