/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_COMPARESTATES_H_
#define TEST_UTIL_COMPARESTATES_H_

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/dot_product.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/IgnoreOtherParameters.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

/// Configuration of the comparestates test.
template <typename MODEL>
class CompareStatesTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CompareStatesTestParameters, Parameters)

 public:
  /// Absolute tolerance of increment norm.
  oops::RequiredParameter<double> tolerance{"tolerance", this};

  /// Configuration of the first state loaded from a file.
  oops::RequiredParameter<eckit::LocalConfiguration> stateConf1{"state1", this};
  /// Configuration of the second state loaded from a file.
  oops::RequiredParameter<eckit::LocalConfiguration> stateConf2{"state2", this};
};

// -----------------------------------------------------------------------------

/// Top-level test parameters.
template <typename MODEL>
class TestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TestParameters, Parameters)

  typedef oops::Geometry<MODEL>              Geometry_;
  typedef typename Geometry_::Parameters_    GeometryParameters_;
  typedef CompareStatesTestParameters<MODEL> CompareStatesTestParameters_;

 public:
  oops::RequiredParameter<CompareStatesTestParameters_> \
    comparestatesTest{"comparestates test", this};
  // Geometry for state1
  oops::RequiredParameter<GeometryParameters_> geometry1{"geometry1", this};
  // Geometry for state2
  oops::RequiredParameter<GeometryParameters_> geometry2{"geometry2", this};
  // The YAML file may also contain options used by other tests; don't treat them as errors.
  oops::IgnoreOtherParameters ignore{this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class CompareStatesFixture : private boost::noncopyable {
 public:
  typedef oops::Geometry<MODEL>              Geometry_;
  typedef CompareStatesTestParameters<MODEL> CompareStatesTestParameters_;

  static const CompareStatesTestParameters_ & test()   {return *getInstance().test_;}
  static const Geometry_                    & resol1() {return *getInstance().resol1_;}
  static const Geometry_                    & resol2() {return *getInstance().resol2_;}
  static void reset() {
    getInstance().test_.reset();
    getInstance().resol1_.reset();
    getInstance().resol2_.reset();
  }

 private:
  static CompareStatesFixture<MODEL>& getInstance() {
    static CompareStatesFixture<MODEL> theCompareStatesFixture;
    return theCompareStatesFixture;
  }

  CompareStatesFixture<MODEL>() {
    TestParameters<MODEL> parameters;
    parameters.validateAndDeserialize(TestEnvironment::config());

    test_   = std::make_unique<CompareStatesTestParameters_>(parameters.comparestatesTest);
    resol1_ = std::make_unique<Geometry_>(parameters.geometry1,
                                         oops::mpi::world(), oops::mpi::myself());
    resol2_ = std::make_unique<Geometry_>(parameters.geometry2,
                                         oops::mpi::world(), oops::mpi::myself());
  }

  ~CompareStatesFixture<MODEL>() {}

  std::unique_ptr<CompareStatesTestParameters_> test_;
  std::unique_ptr<Geometry_>            resol1_;
  std::unique_ptr<Geometry_>            resol2_;
};

// -----------------------------------------------------------------------------
/*! \brief Tests State::geometry() and Geometry copy constructors
 */

template <typename MODEL> void testCompareStates() {
  typedef CompareStatesFixture<MODEL> Test_;
  typedef oops::Geometry<MODEL>       Geometry_;
  typedef oops::State<MODEL>          State_;
  typedef oops::Increment<MODEL>      Increment_;

  const double tol = Test_::test().tolerance;

  State_ xx1(Test_::resol1(), Test_::test().stateConf1);
  State_ xx2(Test_::resol2(), Test_::test().stateConf2);

  // get geometry from xx1
  const Geometry_ & geometry = xx1.geometry();

  //  Create increment (geometry is the same as xx1)
  Increment_ dx(geometry, xx1.variables(), xx1.validTime());
  dx.diff(xx1, xx2);

  EXPECT(oops::is_close_absolute(dx.norm(), 0., tol));
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class CompareStates : public oops::Test {
 public:
  virtual ~CompareStates() {CompareStatesFixture<MODEL>::reset();}

 private:
  std::string testid() const override {return "test::CompareStates<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/CompareStates/testCompareStates")
      { testCompareStates<MODEL>(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_UTIL_COMPARESTATES_H_
