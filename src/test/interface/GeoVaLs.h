/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_GEOVALS_H_
#define TEST_INTERFACE_GEOVALS_H_

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/testing/Test.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/runs/Test.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/IgnoreOtherParameters.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

/// \brief Parameters controlling some of the GeoVaLs test options
template <typename OBS>
class GeoVaLsTestOptionsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeoVaLsTestOptionsParameters, Parameters)
 public:
  oops::RequiredParameter<oops::Variables> vars{"state variables",
                         "variables for GeoVaLs", this};
  oops::RequiredParameter<double> refnorm{"norm",
                         "reference norm for the GeoVaLs test", this};
};

/// \brief Parameters loaded from the input YAML file and used by GeoVaLs test
template <typename OBS>
class GeoVaLsTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeoVaLsTestParameters, Parameters)

 public:
  typedef typename oops::GeoVaLs<OBS>::Parameters_ GeoVaLsParameters_;

  oops::RequiredParameter<GeoVaLsParameters_> geovals{"geovals",
                         "Group of parameters controlling the tested geovals", this};
  oops::RequiredParameter<GeoVaLsTestOptionsParameters<OBS>> geovalstest{"geovals test",
                         "Group of parameters controlling test options", this};
  /// \brief Don't treat the presence of other parameter groups as an error (this makes it
  /// possible to reuse a single YAML file in tests of implementations of multiple oops interfaces).
  oops::IgnoreOtherParameters ignoreOthers{this};
};


// -----------------------------------------------------------------------------
/// \brief Tests test-constructor and print method
template <typename OBS> void testConstructor() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::GeoVaLs<OBS>    GeoVaLs_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    GeoVaLsTestParameters<OBS> params;
    params.validateAndDeserialize(Test_::config(jj));
    std::unique_ptr<GeoVaLs_> geovals(new GeoVaLs_(params.geovals, Test_::obspace()[jj],
                                                   params.geovalstest.value().vars));
    EXPECT(geovals.get());
    oops::Log::info() << "Testing GeoVaLs: " << *geovals << std::endl;
    geovals.reset();
    EXPECT(!geovals.get());
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testUtils() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::GeoVaLs<OBS>    GeoVaLs_;

  const double tol = 1e-6;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    GeoVaLsTestParameters<OBS> params;
    params.validateAndDeserialize(Test_::config(jj));
    GeoVaLs_ gval(params.geovals, Test_::obspace()[jj],
                  params.geovalstest.value().vars);

    const double zz = dot_product(gval, gval);

    oops::Log::trace() << "Testing copy constructor (=) " << std::endl;

    GeoVaLs_ gval2 = gval;

    double zz0 = dot_product(gval2, gval2);

    EXPECT(zz0 == zz);

    oops::Log::trace() << "Testing *= double" << std::endl;

    gval2 *= 2.0;

    const double zz1 = dot_product(gval2, gval2);

    EXPECT(zz1/zz - 4.0 < tol);

    oops::Log::trace() << "Testing += GeoVals" << std::endl;

    gval2 += gval;

    const double zz2 = dot_product(gval2, gval2);

    EXPECT(zz2/zz - 9.0  < tol);

    oops::Log::trace() << "Testing -= GeoVals" << std::endl;

    gval2 -= gval;

    const double zz3 = dot_product(gval2, gval2);

    EXPECT(zz3/zz - 4.0  < tol);

    oops::Log::trace() << "Testing *= GeoVals" << std::endl;

    gval2 = gval;

    gval2 *= gval;

    GeoVaLs_ gval3 = gval2;

    gval3 *= gval;

    const double zz4 = dot_product(gval2, gval2);

    const double zz5 = dot_product(gval3, gval);

    EXPECT(zz4/zz5 - 1.0  < tol);

    oops::Log::trace() << "Testing random" << std::endl;

    gval.random();
    const double zz6 = dot_product(gval, gval);
    EXPECT(zz6 > 0.0);

    oops::Log::trace() << "Testing zero" << std::endl;

    gval.zero();
    const double zz7 = dot_product(gval, gval);
    EXPECT(zz7 == 0.0);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testRead() {
  typedef ObsTestsFixture<OBS> Test_;
  typedef oops::GeoVaLs<OBS>   GeoVaLs_;

  const double tol = 1.0e-9;
  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    GeoVaLsTestParameters<OBS> params;
    params.validateAndDeserialize(Test_::config(jj));
    GeoVaLs_ gval(params.geovals, Test_::obspace()[jj], params.geovalstest.value().vars);

    const double xx = params.geovalstest.value().refnorm;
    const double zz = sqrt(dot_product(gval, gval));

    oops::Log::debug() << "xx: " << std::fixed << std::setprecision(8) << xx << std::endl;
    oops::Log::debug() << "zz: " << std::fixed << std::setprecision(8) << zz << std::endl;

    EXPECT(oops::is_close(xx, zz, tol));
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
class GeoVaLs : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;
 public:
  GeoVaLs() {}
  virtual ~GeoVaLs() {}
 private:
  std::string testid() const override {return "test::GeoVaLs<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/GeoVaLs/testConstructor")
      { testConstructor<OBS>(); });
    ts.emplace_back(CASE("interface/GeoVaLs/testUtils")
      { testUtils<OBS>(); });
    ts.emplace_back(CASE("interface/GeoVaLs/testRead")
      { testRead<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_GEOVALS_H_
