/*
 * (C) Copyright 2020-2022 UCAR
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
#include <utility>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/base/AnalyticInit.h"
#include "oops/base/Geometry.h"
#include "oops/base/GetValues.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/SampledLocations.h"
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/TimeWindow.h"
#include "test/TestEnvironment.h"

namespace test {

// =================================================================================================

template <typename MODEL, typename OBS> class GetValuesFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL> Geometry_;
  typedef oops::SampledLocations<OBS> SampledLocations_;
  typedef oops::Locations<OBS> Locations_;
  typedef oops::State<MODEL>    State_;
  typedef oops::Variables       Variables_;
  typedef util::DateTime        DateTime_;
  typedef util::TimeWindow      TimeWindow_;

 public:
  static const DateTime_         & time()             {return *getInstance().time_;}
  static const TimeWindow_       & timeWindow()       {return *getInstance().timeWindow_;}
  static const Geometry_         & resol()            {return *getInstance().resol_;}
  static const Variables_        & variables()        {return *getInstance().variables_;}
  static const std::vector<size_t> & varsizes()       {return getInstance().varsizes_;}
  static const SampledLocations_ & sampledLocations() {
    // In this test only one location sampling method is used.
    return locations().samplingMethod(0);
  }
  static const Locations_ & locations() {
    return *getInstance().locations_;
  }

  static void reset() {
    getInstance().time_.reset();
    getInstance().timeWindow_.reset();
    getInstance().locations_.reset();
    getInstance().state_.reset();
    getInstance().variables_.reset();
    getInstance().resol_.reset();
  }

 private:
  static GetValuesFixture<MODEL, OBS>& getInstance() {
    static GetValuesFixture<MODEL, OBS> theGetValuesFixture;
    return theGetValuesFixture;
  }

  GetValuesFixture<MODEL, OBS>() {
    // Geometry
    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::world()));

    // Variables
    variables_.reset(new Variables_(TestEnvironment::config(), "variables"));
    varsizes_ = resol_->variableSizes(*variables_);

    // Locations
    const eckit::LocalConfiguration locationsConfig(TestEnvironment::config(), "locations");
    SampledLocations_ sampledLocations(locationsConfig, oops::mpi::world());
    locations_.reset(new Locations_(std::move(sampledLocations)));

    // Assimilation window
    timeWindow_.reset(new TimeWindow_(locationsConfig.getSubConfiguration("time window")));

    // State
    const eckit::LocalConfiguration stateConfig(TestEnvironment::config(), "state");
    state_.reset(new State_(*resol_, stateConfig));
    time_.reset(new DateTime_(state_->validTime()));
  }

  ~GetValuesFixture<MODEL, OBS>() {}

  std::unique_ptr<const DateTime_>        time_;
  std::unique_ptr<const TimeWindow_>      timeWindow_;
  std::unique_ptr<const Geometry_>        resol_;
  std::unique_ptr<const Locations_>       locations_;
  std::unique_ptr<const State_>           state_;
  std::unique_ptr<const Variables_>       variables_;
  std::vector<size_t>                     varsizes_;
};

// -------------------------------------------------------------------------------------------------

/// \brief Test constructor
template <typename MODEL, typename OBS> void testGetValuesConstructor() {
  typedef GetValuesFixture<MODEL, OBS>  Test_;
  typedef oops::GetValues<MODEL, OBS>   GetValues_;

  std::unique_ptr<const GetValues_> getvalues(
      new GetValues_(
          TestEnvironment::config(),
          Test_::resol(),
          Test_::timeWindow(),
          Test_::sampledLocations(),
          Test_::variables()));
  EXPECT(getvalues.get());

  getvalues.reset();
  EXPECT(!getvalues.get());
}

// -------------------------------------------------------------------------------------------------

/// \brief Test GetValues interpolations
///
/// This test constructs two sets of GeoVaLs that should be the same:
/// 1. One GeoVaLs is filled by calling GetValues on a model State initialized from
///    an idealized analytic formula. This interpolates from the analytic function
///    values defined at the model grid points to a set of random locations.
/// 2. The other GeoVaLs is directly filled by evaluating the same analytic formula
///    at the same random locations. This is the expected result if the interpolation in
///    GetValues is of high accuracy.
///
/// These two GeoVaLs are then checked to match within a specified tolerance.
///
/// \note This test requires that the model State has a constructor accepting an "analytic init"
///       config, and a matching AnalyticInit class to fill in the expected GeoVaLs.
template <typename MODEL, typename OBS> void testGetValuesInterpolation() {
  typedef GetValuesFixture<MODEL, OBS>    Test_;
  typedef oops::AnalyticInit<OBS>         AnalyticInit_;
  typedef oops::VariableChange<MODEL>     VariableChange_;
  typedef typename VariableChange_::Parameters_ VariableChangeParameters_;
  typedef oops::State<MODEL>              State_;
  typedef oops::GeoVaLs<OBS>              GeoVaLs_;
  typedef oops::GetValues<MODEL, OBS>     GetValues_;

  const eckit::LocalConfiguration confgen(TestEnvironment::config(), "state");
  const State_ xx(Test_::resol(), confgen);

  eckit::LocalConfiguration chvarconf;  // empty for now
  VariableChangeParameters_ params;
  params.validateAndDeserialize(chvarconf);
  VariableChange_ chvar(params, Test_::resol());
  State_ zz(xx);
  chvar.changeVar(zz, Test_::variables());

  GeoVaLs_ gval(Test_::locations(), Test_::variables(), Test_::varsizes());

  EXPECT(zz.norm() > 0.0);

  // Fill GeoVaLs by calling GetValues to interpolate from the state
  GetValues_ getvalues(TestEnvironment::config(),
                       Test_::resol(),
                       Test_::timeWindow(),
                       Test_::sampledLocations(),
                       Test_::variables());

  oops::PreProcessHelper<MODEL>::preProcessModelData(zz);
  getvalues.initialize(Test_::timeWindow().length());
  getvalues.process(zz);
  getvalues.finalize();
  getvalues.fillGeoVaLs(gval);

  EXPECT(gval.rms() > 0.0);

  // Fill GeoVaLs with exact values
  GeoVaLs_ ref(gval);
  const eckit::LocalConfiguration analyticConf(confgen, "analytic init");
  AnalyticInit_ init(analyticConf);
  init.fillGeoVaLs(Test_::sampledLocations(), ref);

  EXPECT(ref.rms() > 0.0);

  // Compute the difference between the interpolated and exact values
  // and check to see if the errors are within specified tolerance
  const double tol = TestEnvironment::config().getDouble("tolerance interpolation");
  gval -= ref;
  oops::Log::test() << "Normalized rms of the difference: " << gval.normalizedrms(ref) << std::endl;
  EXPECT(gval.normalizedrms(ref) < tol);
}

// -------------------------------------------------------------------------------------------------

/// \brief Test that GetValues on a zero Increment produces a zero GeoVaLs
template <typename MODEL, typename OBS> void testGetValuesTLZeroPert() {
  typedef GetValuesFixture<MODEL, OBS>  Test_;
  typedef oops::GeoVaLs<OBS>            GeoVaLs_;
  typedef oops::GetValues<MODEL, OBS>   GetValues_;
  typedef oops::Increment<MODEL>        Increment_;

  Increment_ dx(Test_::resol(), Test_::variables(), Test_::time());
  dx.zero();
  EXPECT(dx.norm() == 0.0);

  GeoVaLs_ gv(Test_::locations(), Test_::variables(), Test_::varsizes());

  GetValues_ getvalues(TestEnvironment::config(),
                       Test_::resol(),
                       Test_::timeWindow(),
                       Test_::sampledLocations(),
                       Test_::variables(),
                       Test_::variables());  // linear variables


  // Test passing zeros forward
  const util::Duration windowlength = Test_::timeWindow().length();
  oops::PreProcessHelper<MODEL>::preProcessModelData(dx);
  getvalues.initializeTL(windowlength);
  getvalues.processTL(dx);
  getvalues.finalizeTL();
  getvalues.fillGeoVaLsTL(gv);

  EXPECT(dx.norm() == 0.0);
  EXPECT(gv.rms() == 0.0);

  // Test adjoint of passing zeros
  getvalues.fillGeoVaLsAD(gv);
  getvalues.finalizeAD(windowlength);
  getvalues.processAD(dx);
  getvalues.initializeAD();
  oops::PreProcessHelper<MODEL>::preProcessModelDataAD(dx);
  dx.synchronizeFields();

  EXPECT(dx.norm() == 0.0);
  EXPECT(gv.rms() == 0.0);
}

// -------------------------------------------------------------------------------------------------

/// \brief Test that GetValues is a linear function of the input Increment
template <typename MODEL, typename OBS> void testGetValuesLinearity() {
  typedef GetValuesFixture<MODEL, OBS>  Test_;
  typedef oops::GeoVaLs<OBS>            GeoVaLs_;
  typedef oops::GetValues<MODEL, OBS>   GetValues_;
  typedef oops::Increment<MODEL>        Increment_;

  const double zz = 3.1415;  // arbitrary scalar factor

  Increment_ dx1(Test_::resol(), Test_::variables(), Test_::time());
  dx1.random();
  Increment_ dx2(dx1);

  EXPECT(dx1.norm() > 0.0);
  EXPECT(dx2.norm() > 0.0);

  GeoVaLs_ gv1(Test_::locations(), Test_::variables(), Test_::varsizes());
  GeoVaLs_ gv2(Test_::locations(), Test_::variables(), Test_::varsizes());

  const util::Duration windowlength = Test_::timeWindow().length();
  GetValues_ getvalues(TestEnvironment::config(),
                       Test_::resol(),
                       Test_::timeWindow(),
                       Test_::sampledLocations(),
                       Test_::variables(),
                       Test_::variables());  // linear variables

  // Compute geovals
  oops::PreProcessHelper<MODEL>::preProcessModelData(dx1);
  getvalues.initializeTL(windowlength);
  getvalues.processTL(dx1);
  getvalues.finalizeTL();
  getvalues.fillGeoVaLsTL(gv1);

  gv1 *= zz;
  dx2 *= zz;

  // Compute geovals
  oops::PreProcessHelper<MODEL>::preProcessModelData(dx2);
  getvalues.initializeTL(windowlength);
  getvalues.processTL(dx2);
  getvalues.finalizeTL();
  getvalues.fillGeoVaLsTL(gv2);

  const double tol = TestEnvironment::config().getDouble("tolerance linearity", 1.0e-11);
  GeoVaLs_ gv_diff = gv1;
  gv_diff -= gv2;
  EXPECT(gv_diff.rms() < tol);
}

// -------------------------------------------------------------------------------------------------

/// \brief Test the adjoint of GetValues
template <typename MODEL, typename OBS> void testGetValuesAdjoint() {
  typedef GetValuesFixture<MODEL, OBS>  Test_;
  typedef oops::GeoVaLs<OBS>            GeoVaLs_;
  typedef oops::GetValues<MODEL, OBS>   GetValues_;
  typedef oops::Increment<MODEL>        Increment_;

  Increment_ dx_in(Test_::resol(), Test_::variables(), Test_::time());
  Increment_ dx_out(Test_::resol(), Test_::variables(), Test_::time());

  GeoVaLs_ gv_out(Test_::locations(), Test_::variables(), Test_::varsizes());

  GetValues_ getvalues(TestEnvironment::config(),
                       Test_::resol(),
                       Test_::timeWindow(),
                       Test_::sampledLocations(),
                       Test_::variables(),
                       Test_::variables());  // linear variables

  // Tangent linear
  const util::Duration windowlength = Test_::timeWindow().length();
  dx_in.random();
  EXPECT(dx_in.norm() > 0.0);
  oops::PreProcessHelper<MODEL>::preProcessModelData(dx_in);
  getvalues.initializeTL(windowlength);
  getvalues.processTL(dx_in);
  getvalues.finalizeTL();
  getvalues.fillGeoVaLsTL(gv_out);
  EXPECT(gv_out.rms() > 0.0);

  // Adjoint
  GeoVaLs_ gv_in(gv_out);
  gv_in.random();
  EXPECT(gv_in.rms() > 0.0);
  dx_out.zero();
  getvalues.fillGeoVaLsAD(gv_in);
  getvalues.finalizeAD(windowlength);
  getvalues.processAD(dx_out);
  getvalues.initializeAD();
  oops::PreProcessHelper<MODEL>::preProcessModelDataAD(dx_out);
  dx_out.synchronizeFields();
  EXPECT(dx_out.norm() > 0.0);

  // Dot products
  const double dot1 = dot_product(dx_in, dx_out);
  const double dot2 = dot_product(gv_in, gv_out);
  const double tol = TestEnvironment::config().getDouble("tolerance AD", 1.0e-11);
  EXPECT(oops::is_close(dot1, dot2, tol));

  oops::Log::test() << "Dot Product <dx, M^Tgv> = " << dot1 << std::endl;
  oops::Log::test() << "Dot Product <gv, M  dx> = " << dot2 << std::endl;
  oops::Log::test() << "Relative diff: " << (dot1-dot2)/dot1 << std::endl;
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL, typename OBS>
class GetValues : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~GetValues() {GetValuesFixture<MODEL, OBS>::reset();}

 private:
  std::string testid() const override {return "test::GetValues<" + MODEL::name() +
                                              ", " + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/GetValues/testGetValuesConstructor")
      { testGetValuesConstructor<MODEL, OBS>(); });
    ts.emplace_back(CASE("interface/GetValues/testGetValuesInterpolation")
      { testGetValuesInterpolation<MODEL, OBS>(); });
    ts.emplace_back(CASE("interface/GetValues/testGetValuesTLZeroPert")
      { testGetValuesTLZeroPert<MODEL, OBS>(); });
    ts.emplace_back(CASE("interface/GetValues/testGetValuesLinearity")
      { testGetValuesLinearity<MODEL, OBS>(); });
    ts.emplace_back(CASE("interface/GetValues/testGetValuesAdjoint")
      { testGetValuesAdjoint<MODEL, OBS>(); });
  }

  void clear() const override {}
};

// =================================================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GETVALUES_H_
