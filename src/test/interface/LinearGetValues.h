/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_LINEARGETVALUES_H_
#define TEST_INTERFACE_LINEARGETVALUES_H_

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
#include "oops/interface/Geometry.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/GetValues.h"
#include "oops/interface/Increment.h"
#include "oops/interface/LinearGetValues.h"
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

template <typename MODEL, typename OBS> class LinearGetValuesFixture : private boost::noncopyable {
  typedef eckit::LocalConfiguration          LocalConfig_;
  typedef oops::GeoVaLs<OBS>                 GeoVaLs_;
  typedef oops::Geometry<MODEL>              Geometry_;
  typedef oops::GetValues<MODEL, OBS>        GetValues_;
  typedef oops::LinearGetValues<MODEL, OBS>  LinearGetValues_;
  typedef oops::Locations<OBS>               Locations_;
  typedef oops::State<MODEL>                 State_;
  typedef oops::Variables                    Variables_;
  typedef util::DateTime                     DateTime_;

 public:
  static const DateTime_         & time()            {return *getInstance().time_;}
  static const DateTime_         & timebeg()         {return *getInstance().timebeg_;}
  static const DateTime_         & timeend()         {return *getInstance().timeend_;}
  static const GeoVaLs_          & geovals()         {return *getInstance().geovals_;}
  static const Geometry_         & resol()           {return *getInstance().resol_;}
  static const GetValues_        & getvalues()       {return *getInstance().getvalues_;}
  static const LinearGetValues_  & lineargetvalues() {return *getInstance().lineargetvalues_;}
  static const LocalConfig_      & testconf()        {return *getInstance().testconf_;}
  static const Locations_        & locs()            {return *getInstance().locs_;}
  static const State_            & state()           {return *getInstance().state_;}
  static const Variables_        & statevars()       {return *getInstance().statevars_;}
  static const Variables_        & geovalvars()      {return *getInstance().geovalvars_;}
  static const std::vector<size_t> & geovalvarsizes() {return getInstance().geovalvarsizes_;}

  static void reset() {
    getInstance().lineargetvalues_.reset();
    getInstance().time_.reset();
    getInstance().statevars_.reset();
    getInstance().state_.reset();
    getInstance().getvalues_.reset();
    getInstance().geovals_.reset();
    getInstance().timeend_.reset();
    getInstance().timebeg_.reset();
    getInstance().locs_.reset();
    getInstance().geovalvars_.reset();
    getInstance().resol_.reset();
    getInstance().testconf_.reset();
  }

 private:
  static LinearGetValuesFixture<MODEL, OBS>& getInstance() {
    static LinearGetValuesFixture<MODEL, OBS> theLinearGetValuesFixture;
    return theLinearGetValuesFixture;
  }

  LinearGetValuesFixture<MODEL, OBS>() {
    testconf_.reset(new LocalConfig_(TestEnvironment::config(), "linear getvalues test"));

    // Geometry
    const LocalConfig_ resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::world()));

    // Variables
    geovalvars_.reset(new Variables_(TestEnvironment::config(), "state variables"));
    geovalvarsizes_ = resol_->variableSizes(*geovalvars_);

    // Locations
    const LocalConfig_ locsConfig(TestEnvironment::config(), "locations");
    locs_.reset(new Locations_(locsConfig, oops::mpi::world()));

    // Window times
    timebeg_.reset(new DateTime_(locsConfig.getString("window begin")));
    timeend_.reset(new DateTime_(locsConfig.getString("window end")));

    // GeoVaLs
    geovals_.reset(new GeoVaLs_(*locs_, *geovalvars_, geovalvarsizes_));

    // Nonlinear GetValues
    LocalConfig_ getvaluesConfig;
    if (TestEnvironment::config().has("get values"))
      getvaluesConfig = eckit::LocalConfiguration(TestEnvironment::config(), "get values");
    getvalues_.reset(new GetValues_( *resol_, *locs_, getvaluesConfig));

    // State
    const LocalConfig_ stateConfig(TestEnvironment::config(), "background");
    state_.reset(new State_(*resol_, stateConfig));
    statevars_.reset(new Variables_(state_->variables()));

    // Valid time
    time_.reset(new DateTime_(state_->validTime()));

    // LinearGetValues
    LocalConfig_ linearGetValuesConfig;
    if (TestEnvironment::config().has("linear get values"))
      linearGetValuesConfig =
        eckit::LocalConfiguration(TestEnvironment::config(), "linear get values");
    lineargetvalues_.reset(new LinearGetValues_(*resol_, *locs_,
                                                linearGetValuesConfig));

    // Set trajectory
    GeoVaLs_ gvtraj(*locs_, *geovalvars_, geovalvarsizes_);
    lineargetvalues_->setTrajectory(*state_, *timebeg_, *timeend_, gvtraj);
  }

  ~LinearGetValuesFixture<MODEL, OBS>() {}

  std::unique_ptr<const DateTime_>        time_;
  std::unique_ptr<const DateTime_>        timebeg_;
  std::unique_ptr<const DateTime_>        timeend_;
  std::unique_ptr<const GeoVaLs_>         geovals_;
  std::unique_ptr<const Geometry_>        resol_;
  std::unique_ptr<const GetValues_>       getvalues_;
  std::unique_ptr<LinearGetValues_>       lineargetvalues_;
  std::unique_ptr<const LocalConfig_>     testconf_;
  std::unique_ptr<const Locations_>       locs_;
  std::unique_ptr<const State_>           state_;
  std::unique_ptr<const Variables_>       statevars_;
  std::unique_ptr<const Variables_>       geovalvars_;
  std::vector<size_t>                     geovalvarsizes_;
};

// =================================================================================================
/// \brief tests constructor and print method
template <typename MODEL, typename OBS> void testLinearGetValuesConstructor() {
  typedef LinearGetValuesFixture<MODEL, OBS>  Test_;
  typedef oops::LinearGetValues<MODEL, OBS>   LinearGetValues_;

  eckit::LocalConfiguration linearGetValuesConfig;
  if (TestEnvironment::config().has("linear get values"))
    linearGetValuesConfig =
      eckit::LocalConfiguration(TestEnvironment::config(), "linear get values");

  std::unique_ptr<const LinearGetValues_>
    lineargetvalues(new LinearGetValues_(Test_::resol(),
                                         Test_::locs(),
                                         linearGetValuesConfig));
  EXPECT(lineargetvalues.get());
  oops::Log::test() << "Testing LinearGetValues: " << *lineargetvalues << std::endl;
  lineargetvalues.reset();
  EXPECT(!lineargetvalues.get());
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL, typename OBS> void testLinearGetValuesZeroPert() {
  typedef LinearGetValuesFixture<MODEL, OBS>  Test_;
  typedef oops::GeoVaLs<OBS>                  GeoVaLs_;
  typedef oops::Increment<MODEL>              Increment_;

  Increment_ dx(Test_::resol(), Test_::statevars(), Test_::time());
  dx.zero();

  GeoVaLs_ gv(Test_::locs(), Test_::geovalvars(), Test_::geovalvarsizes());

  EXPECT(dx.norm() == 0.0);

  Test_::lineargetvalues().fillGeoVaLsTL(dx, Test_::timebeg(), Test_::timeend(), gv);

  EXPECT(dx.norm() == 0.0);
  EXPECT(gv.rms() == 0.0);

  Test_::lineargetvalues().fillGeoVaLsAD(dx, Test_::timebeg(), Test_::timeend(), gv);

  EXPECT(dx.norm() == 0.0);
  EXPECT(gv.rms() == 0.0);
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL, typename OBS> void testLinearGetValuesLinearity() {
  typedef LinearGetValuesFixture<MODEL, OBS>  Test_;
  typedef oops::GeoVaLs<OBS>                  GeoVaLs_;
  typedef oops::Increment<MODEL>              Increment_;

  const double zz = 3.1415;

  Increment_ dx1(Test_::resol(), Test_::statevars(), Test_::time());
  dx1.random();
  Increment_ dx2(dx1);

  EXPECT(dx1.norm() > 0.0);
  EXPECT(dx2.norm() > 0.0);

  GeoVaLs_ gv1(Test_::locs(), Test_::geovalvars(), Test_::geovalvarsizes());
  GeoVaLs_ gv2(Test_::locs(), Test_::geovalvars(), Test_::geovalvarsizes());

  // Compute geovals
  Test_::lineargetvalues().fillGeoVaLsTL(dx1, Test_::timebeg(), Test_::timeend(), gv1);

  gv1 *= zz;
  dx2 *= zz;

  // Compute geovals
  Test_::lineargetvalues().fillGeoVaLsTL(dx2, Test_::timebeg(), Test_::timeend(), gv2);

  const double tol = Test_::testconf().getDouble("tolerance linearity", 1.0e-11);
  EXPECT(oops::is_close(gv1.rms(), gv2.rms(), tol));
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL, typename OBS> void testLinearGetValuesLinearApproximation() {
  typedef LinearGetValuesFixture<MODEL, OBS>  Test_;
  typedef oops::GeoVaLs<OBS>                  GeoVaLs_;
  typedef oops::Increment<MODEL>              Increment_;
  typedef oops::State<MODEL>                  State_;

  const unsigned int ntest = Test_::testconf().getInt("iterations TL", 10);
  double zz = Test_::testconf().getDouble("first multiplier TL", 1.0e-2);

  // Compute nonlinear geovals
  State_ xx0(Test_::state());
  GeoVaLs_ gv0(Test_::locs(), Test_::geovalvars(), Test_::geovalvarsizes());
  Test_::getvalues().fillGeoVaLs(xx0, Test_::timebeg(), Test_::timeend(), gv0);

  // Run tangent linear
  Increment_ dx0(Test_::resol(), Test_::statevars(), Test_::time());
  dx0.random();

  std::vector<double> errors;
  for (unsigned int jtest = 0; jtest < ntest; ++jtest) {
    Increment_ dxx(dx0);
    State_ xx(xx0);
    GeoVaLs_  gv(Test_::locs(), Test_::geovalvars(), Test_::geovalvarsizes());
    GeoVaLs_ dgv(Test_::locs(), Test_::geovalvars(), Test_::geovalvarsizes());
    dxx *= zz;
    xx += dxx;

    // Nonlinear
    Test_::getvalues().fillGeoVaLs(xx, Test_::timebeg(), Test_::timeend(), gv);

    // Tangent linear
    Test_::lineargetvalues().fillGeoVaLsTL(dxx, Test_::timebeg(), Test_::timeend(), dgv);

    GeoVaLs_ gvdiff(gv);
    gvdiff -= gv0;

    // Print the norms
    const double nlnorm = gvdiff.rms();
    const double tlnorm = dgv.rms();

    // Subtract the tlm pert
    gvdiff -= dgv;

    double testdot = dot_product(gvdiff, gvdiff);
    errors.push_back(testdot);

    oops::Log::test() << "  ||g(x+dx) - g(x)|| = " << std::setprecision(16) << nlnorm << std::endl;
    oops::Log::test() << "             ||Gdx|| = " << std::setprecision(16) << tlnorm << std::endl;
    oops::Log::test() << "||g(x+dx)-g(x)-Gdx|| = " << std::setprecision(16) << testdot << std::endl;

    zz /= 10.0;
  }

  // Analyze results
  const double approx = *std::min_element(errors.begin(), errors.end());
  oops::Log::test() << "Test LinearGetValuesTL min error = " << approx << std::endl;
  const double tol = Test_::testconf().getDouble("tolerance TL", 1.0e-11);
  EXPECT(approx < tol);
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL, typename OBS> void testLinearGetValuesAdjoint() {
  typedef LinearGetValuesFixture<MODEL, OBS>  Test_;
  typedef oops::GeoVaLs<OBS>                  GeoVaLs_;
  typedef oops::Increment<MODEL>              Increment_;

  Increment_ dx_in(Test_::resol(), Test_::statevars(), Test_::time());
  Increment_ dx_ou(Test_::resol(), Test_::statevars(), Test_::time());

  GeoVaLs_ gv_ou(Test_::locs(), Test_::geovalvars(), Test_::geovalvarsizes());

  // Tangent linear
  dx_in.random();
  EXPECT(dx_in.norm() > 0.0);
  Test_::lineargetvalues().fillGeoVaLsTL(dx_in, Test_::timebeg(), Test_::timeend(), gv_ou);
  EXPECT(gv_ou.rms() > 0.0);

  // Adjoint
  GeoVaLs_ gv_in(gv_ou);
  gv_in.random();  // No order dependency but need a copy
  EXPECT(gv_in.rms() > 0.0);
  dx_ou.zero();
  Test_::lineargetvalues().fillGeoVaLsAD(dx_ou, Test_::timebeg(), Test_::timeend(), gv_in);
  EXPECT(dx_ou.norm() > 0.0);

  // Dot products
  const double dot1 = dot_product(dx_in, dx_ou);
  const double dot2 = dot_product(gv_in, gv_ou);
  const double tol = Test_::testconf().getDouble("tolerance AD", 1.0e-11);
  EXPECT(oops::is_close(dot1, dot2, tol));

  oops::Log::test() << "Dot Product <dx, M^Tgv> = " << std::setprecision(16) << dot1 << std::endl;
  oops::Log::test() << "Dot Product <gv, M  dx> = " << std::setprecision(16) << dot2 << std::endl;
  oops::Log::test() << "Relative diff: " << std::setprecision(16) << (dot1-dot2)/dot1 << std::endl;
}

// =================================================================================================

template <typename MODEL, typename OBS>
class LinearGetValues : public oops::Test {
 public:
  LinearGetValues() {}
  virtual ~LinearGetValues() {LinearGetValuesFixture<MODEL, OBS>::reset();}

 private:
  std::string testid() const override {return "test::LinearGetValues<" + MODEL::name() + ", "
                                                              + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/LinearGetValues/testLinearGetValuesConstructor")
      { testLinearGetValuesConstructor<MODEL, OBS>(); });
    ts.emplace_back(CASE("interface/LinearGetValues/testLinearGetValuesZeroPert")
      { testLinearGetValuesZeroPert<MODEL, OBS>(); });
    ts.emplace_back(CASE("interface/LinearGetValues/testLinearGetValuesLinearity")
      { testLinearGetValuesLinearity<MODEL, OBS>(); });
    ts.emplace_back(CASE("interface/LinearGetValues/testLinearGetValuesLinearApproximation")
      { testLinearGetValuesLinearApproximation<MODEL, OBS>(); });
    ts.emplace_back(CASE("interface/LinearGetValues/testLinearGetValuesAdjoint")
      { testLinearGetValuesAdjoint<MODEL, OBS>(); });
  }

  void clear() const override {}
};

// =================================================================================================

}  // namespace test

#endif  // TEST_INTERFACE_LINEARGETVALUES_H_
