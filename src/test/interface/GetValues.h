/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_GETVALUES_H_
#define TEST_INTERFACE_GETVALUES_H_

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Variables.h"
#include "oops/interface/AnalyticInit.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/GetValues.h"
#include "oops/interface/Locations.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "test/TestEnvironment.h"

namespace test {

// =================================================================================================

template <typename MODEL, typename OBS> class GetValuesFixture : private boost::noncopyable {
  typedef eckit::LocalConfiguration    LocalConfig_;
  typedef oops::GeoVaLs<OBS>           GeoVaLs_;
  typedef oops::Geometry<MODEL>        Geometry_;
  typedef oops::GetValues<MODEL, OBS>  GetValues_;
  typedef oops::Locations<OBS>         Locations_;
  typedef oops::Variables              Variables_;
  typedef util::DateTime               DateTime_;

 public:
  static const DateTime_         & timebeg()         {return *getInstance().timebeg_;}
  static const DateTime_         & timeend()         {return *getInstance().timeend_;}
  static const GeoVaLs_          & geovals()         {return *getInstance().geovals_;}
  static const Geometry_         & resol()           {return *getInstance().resol_;}
  static const GetValues_        & getvalues()       {return *getInstance().getvalues_;}
  static const LocalConfig_      & testconf()        {return *getInstance().testconf_;}
  static const Locations_        & locs()            {return *getInstance().locs_;}
  static const Variables_        & geovalvars()      {return *getInstance().geovalvars_;}

 private:
  static GetValuesFixture<MODEL, OBS>& getInstance() {
    static GetValuesFixture<MODEL, OBS> theGetValuesFixture;
    return theGetValuesFixture;
  }

  GetValuesFixture<MODEL, OBS>() {
    testconf_.reset(new LocalConfig_(TestEnvironment::config(), "getvalues test"));

    // Geometry
    const LocalConfig_ resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::world()));

    // Variables
    geovalvars_.reset(new Variables_(TestEnvironment::config(), "state variables"));

    // Locations
    const LocalConfig_ locsConfig(TestEnvironment::config(), "locations");
    locs_.reset(new Locations_(locsConfig, oops::mpi::world()));

    // Window times
    timebeg_.reset(new DateTime_(locsConfig.getString("window begin")));
    timeend_.reset(new DateTime_(locsConfig.getString("window end")));

    // GeoVaLs
    geovals_.reset(new GeoVaLs_(*locs_, *geovalvars_));

    // GetValues
    getvalues_.reset(new GetValues_(*resol_, *locs_));
  }

  ~GetValuesFixture<MODEL, OBS>() {}

  std::unique_ptr<const DateTime_>        timebeg_;
  std::unique_ptr<const DateTime_>        timeend_;
  std::unique_ptr<const GeoVaLs_>         geovals_;
  std::unique_ptr<const Geometry_>        resol_;
  std::unique_ptr<const GetValues_>       getvalues_;
  std::unique_ptr<const LocalConfig_>     testconf_;
  std::unique_ptr<const Locations_>       locs_;
  std::unique_ptr<const Variables_>       geovalvars_;
};

// =================================================================================================
/// \brief tests constructor and print method
template <typename MODEL, typename OBS> void testGetValuesConstructor() {
  typedef GetValuesFixture<MODEL, OBS>  Test_;
  typedef oops::GetValues<MODEL, OBS>   GetValues_;

  std::unique_ptr<const GetValues_> GetValues(new GetValues_(Test_::resol(), Test_::locs()));
  EXPECT(GetValues.get());
  oops::Log::test() << "Testing GetValues: " << *GetValues << std::endl;
  GetValues.reset();
  EXPECT(!GetValues.get());
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL, typename OBS> void testGetValuesMultiWindow() {
  typedef GetValuesFixture<MODEL, OBS>    Test_;
  typedef oops::GeoVaLs<OBS>              GeoVaLs_;
  typedef oops::State<MODEL>         State_;

  const util::Duration windowlength = Test_::timeend() - Test_::timebeg();
  const util::DateTime timemid = Test_::timebeg() + windowlength/2;

  const eckit::LocalConfiguration confgen(Test_::testconf(), "state generate");
  const State_ xx(Test_::resol(), confgen);

  EXPECT(xx.norm() > 0.0);

  GeoVaLs_ gv1(Test_::locs(), Test_::geovalvars());
  GeoVaLs_ gv2(Test_::locs(), Test_::geovalvars());

  // Compute all geovals together
  Test_::getvalues().fillGeoVaLs(xx, Test_::timebeg(), Test_::timeend(), gv1);

  // Compute all geovals as two subwindows
  Test_::getvalues().fillGeoVaLs(xx, Test_::timebeg(), timemid, gv2);
  Test_::getvalues().fillGeoVaLs(xx, timemid, Test_::timeend(), gv2);

  EXPECT(gv1.rms() > 0.0);
  EXPECT(gv2.rms() > 0.0);
  EXPECT(gv1.rms() == gv2.rms());
}

// -------------------------------------------------------------------------------------------------

/*! \brief Interpolation test
 *
 * \details **testGetValuesInterpolation()** tests the creation of an
 * analytic state for a given model.  The conceptual steps are as follows:
 * 1. Initialize the JEDI State object based on idealized analytic formulae
 * 2. Interpolate the State variables onto selected "observation" locations
 *    using the getValues() method of the State object.  The result is
 *    placed in a JEDI GeoVaLs object
 * 3. Compute the correct solution by applying the analytic formulae directly
 *    at the observation locations.
 * 4. Assess the accuracy of the interpolation by comparing the interpolated
 *    values from Step 2 with the exact values from Step 3
 *
 * The interpolated state values are compared to the analytic solution for
 * a series of **locations** which includes values optionally specified by the
 * user in the "state test" section of the config file in addition to a
 * randomly-generated list of **nrandom** random locations.  nrandom is also
 * specified by the user in the "state test" section of the config file, as is the
 * (nondimensional) tolerence level (**interpolation tolerance**) to be used for the tests.
 *
 * Relevant parameters in the **State* section of the config file include
 *
 * * **norm-gen** Normalization test for the generated State
 * * **interpolation tolerance** tolerance for the interpolation test
 *
 * \date April, 2018: M. Miesch (JCSDA) adapted a preliminary version in the
 * feature/interp branch
 *
 * \warning Since this model compares the interpolated state values to an exact analytic
 * solution, it requires that the "analytic_init" option be implemented in the model and
 * selected in the "state.state generate" section of the config file.
 */

template <typename MODEL, typename OBS> void testGetValuesInterpolation() {
  typedef GetValuesFixture<MODEL, OBS>    Test_;
  typedef oops::AnalyticInit<OBS>         AnalyticInit_;
  typedef oops::State<MODEL>              State_;
  typedef oops::GeoVaLs<OBS>              GeoVaLs_;

  const eckit::LocalConfiguration confgen(Test_::testconf(), "state generate");
  const State_ xx(Test_::resol(), confgen);

  // Interpolation tolerance
  double interp_tol = Test_::testconf().getDouble("interpolation tolerance");

  // Ceate a GeoVaLs object from locs and vars
  GeoVaLs_ gval(Test_::locs(), Test_::geovalvars());

  EXPECT(xx.norm() > 0.0);

  // Execute the interpolation
  Test_::getvalues().fillGeoVaLs(xx, Test_::timebeg(), Test_::timeend(), gval);

  EXPECT(gval.rms() > 0.0);
  oops::Log::debug() << "RMS GeoVaLs: " << gval.rms() << std::endl;

  // Now create another GeoVaLs object that contains the exact analytic solutions.
  GeoVaLs_ ref(gval);

  AnalyticInit_ init(confgen);
  init.fillGeoVaLs(Test_::locs(), ref);

  EXPECT(ref.rms() > 0.0);
  oops::Log::debug() << "RMS reference GeoVaLs: " << ref.rms() << std::endl;

  // Compute the difference between the interpolated and exact values
  gval -= ref;

  // And check to see if the errors are within specified tolerance
  oops::Log::test() << "Normalized rms of the difference: " << gval.normalizedrms(ref) << std::endl;
  EXPECT(gval.normalizedrms(ref) < interp_tol);
}

// =================================================================================================

template <typename MODEL, typename OBS>
class GetValues : public oops::Test {
 public:
  GetValues() {}
  virtual ~GetValues() {}

 private:
  std::string testid() const override {return "test::GetValues<" + MODEL::name() +
                                              ", " + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/GetValues/testGetValuesConstructor")
      { testGetValuesConstructor<MODEL, OBS>(); });
    ts.emplace_back(CASE("interface/GetValues/testGetValuesMultiWindow")
      { testGetValuesMultiWindow<MODEL, OBS>(); });
    ts.emplace_back(CASE("interface/GetValues/testGetValuesInterpolation")
      { testGetValuesInterpolation<MODEL, OBS>(); });
  }

  void clear() const override {}
};

// =================================================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GETVALUES_H_
