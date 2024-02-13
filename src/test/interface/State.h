/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_STATE_H_
#define TEST_INTERFACE_STATE_H_

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
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
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

/// Configuration of the state test.
template <typename MODEL>
class StateTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StateTestParameters, Parameters)

 public:
  /// Relative tolerance of norm comparisons.
  oops::RequiredParameter<double> tolerance{"tolerance", this};
  /// Validity time for states loaded from a file and generated on the fly.
  oops::RequiredParameter<util::DateTime> date{"date", this};

  /// Configuration of the state loaded from a file.
  oops::RequiredParameter<eckit::LocalConfiguration> statefile{"statefile", this};
  /// Expected norm of the state loaded from a file.
  oops::RequiredParameter<double> normFile{"norm file", this};

  /// Configuration of the state generated on the fly.
  oops::OptionalParameter<eckit::LocalConfiguration> stateGenerate{"state generate", this};
  /// Expected norm of the state generated on the fly.
  ///
  /// This option must be present if `state generate` is.
  oops::OptionalParameter<double> normGeneratedState{"norm generated state", this};

  oops::OptionalParameter<eckit::LocalConfiguration> writeReadTest{"write then read test", this};

  /// Flag indicating whether to run the test of the variable change State constructor
  oops::Parameter<bool> testVarConstructor{"test variable change constructor", true, this};

  /// Variables to pass to the variable change State constructor in the test
  oops::OptionalParameter<oops::Variables> toVariables{"construct to variables", this};
};

// -----------------------------------------------------------------------------

/// Top-level test parameters.
template <typename MODEL>
class TestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TestParameters, Parameters)

  typedef oops::Geometry<MODEL>           Geometry_;
  typedef typename Geometry_::Parameters_ GeometryParameters_;
  typedef StateTestParameters<MODEL>      StateTestParameters_;

 public:
  oops::RequiredParameter<StateTestParameters_> stateTest{"state test", this};
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};
  // The YAML file may also contain options used by other tests; don't treat them as errors.
  oops::IgnoreOtherParameters ignore{this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class StateFixture : private boost::noncopyable {
 public:
  typedef oops::Geometry<MODEL>      Geometry_;
  typedef StateTestParameters<MODEL> StateTestParameters_;

  static const StateTestParameters_ & test()  {return *getInstance().test_;}
  static const Geometry_            & resol() {return *getInstance().resol_;}
  static void reset() {
    getInstance().resol_.reset();
    getInstance().test_.reset();
  }

 private:
  static StateFixture<MODEL>& getInstance() {
    static StateFixture<MODEL> theStateFixture;
    return theStateFixture;
  }

  StateFixture<MODEL>() {
    TestParameters<MODEL> parameters;
    parameters.validateAndDeserialize(TestEnvironment::config());

    test_ = std::make_unique<StateTestParameters_>(parameters.stateTest);
    resol_ = std::make_unique<Geometry_>(parameters.geometry,
                                         oops::mpi::world(), oops::mpi::myself());
  }

  ~StateFixture<MODEL>() {}

  std::unique_ptr<StateTestParameters_> test_;
  std::unique_ptr<Geometry_>            resol_;
};

// -----------------------------------------------------------------------------
/// \brief tests constructors and print method
template <typename MODEL> void testStateConstructors() {
  typedef StateFixture<MODEL>   Test_;
  typedef oops::State<MODEL>    State_;

  const double norm = Test_::test().normFile;
  const double tol = Test_::test().tolerance;
  const util::DateTime vt(Test_::test().date);

// Test main constructor
  std::unique_ptr<State_> xx1(new State_(Test_::resol(), Test_::test().statefile));

  EXPECT(xx1.get());
  oops::Log::test() << "Printing State from yaml: " << *xx1 << std::endl;
  const double norm1 = xx1->norm();
  EXPECT(norm1 != 0);
  EXPECT(oops::is_close(norm1, norm, tol));
  EXPECT(xx1->validTime() == vt);

// Test copy constructor
  std::unique_ptr<State_> xx2(new State_(*xx1));
  EXPECT(xx2.get());
  EXPECT(oops::is_close(xx2->norm(), norm, tol));
  EXPECT(xx2->validTime() == vt);

// Destruct copy
  xx2.reset();
  EXPECT(!xx2.get());

// Recompute initial norm to make sure nothing bad happened
  const double norm2 = xx1->norm();
  EXPECT(norm1 == norm2);

// Test State(const Geometry_ &, const Variables &, const util::DateTime &) constructor
  oops::Variables vars(xx1->variables());
  State_ xx3(Test_::resol(), vars, vt);
  oops::Log::test() << "Printing empty State: " << xx3 << std::endl;
  EXPECT(xx3.norm() == 0);
  EXPECT(xx3.validTime() == vt);
  EXPECT(xx3.variables() == vars);

// Test State(const Geometry_ &, const State &) constructor
  State_ xx4(Test_::resol(), *xx1);
  EXPECT(oops::is_close(xx4.norm(), norm, tol));
  EXPECT(xx4.validTime() == vt);
  EXPECT(xx4.variables() == xx1->variables());

// Test State(const Variables &, const State &) constructor (unless told not to)
  const bool testVarChangeConstructor = Test_::test().testVarConstructor.value();
  if (testVarChangeConstructor)
  {
    EXPECT(Test_::test().toVariables.value() != boost::none);
    const auto & toVars = Test_::test().toVariables.value().value();
    EXPECT(toVars.size() > 0);
    State_ xx5(toVars, *xx1);
    EXPECT(xx5.norm() > 0.0);
    EXPECT(xx5.validTime() == vt);
    EXPECT(xx5.variables() == toVars);
  }
}

// -----------------------------------------------------------------------------
/// \brief tests constructors and print method
template <typename MODEL> void testStateAtlasInterface() {
  typedef StateFixture<MODEL>   Test_;
  typedef oops::State<MODEL>    State_;

  const bool testAtlas = TestEnvironment::config().getBool("test atlas interface", true);
  if (!testAtlas) { return; }

  const oops::Geometry<MODEL> & geom = Test_::resol();

  State_ xx(geom, Test_::test().statefile);
  const oops::Variables & vars = xx.variables();

  atlas::FieldSet fset{};
  xx.toFieldSet(fset);

  // Check Fields in FieldSet
  const int nvars = vars.size();
  EXPECT(fset.size() == nvars);
  for (int v = 0; v < nvars; ++v) {
    const atlas::Field & f = fset[v];
    EXPECT(f.valid());
    EXPECT(f.functionspace() == geom.functionSpace());
    EXPECT(!f.dirty());
    EXPECT(f.rank() == 2);
    EXPECT(f.shape(0) == geom.functionSpace().lonlat().shape(0));
    EXPECT(f.datatype() == atlas::array::DataType::create<double>());
  }

  // Check haloExchange is no-op, i.e., halos are up-to-date
  atlas::FieldSet fset2 = util::copyFieldSet(fset);
  for (int v = 0; v < nvars; ++v) {
    fset2[v].set_dirty();
  }
  fset2.haloExchange();
  EXPECT(util::compareFieldSets(fset, fset2));

  // Check fromFieldSet repopulates State in the same way
  const util::DateTime t(Test_::test().date);
  State_ yy(geom, vars, t);
  yy.fromFieldSet(fset);

  // How to check xx and yy are the same state?
  // - They might not be, because States can include fields that aren't part of the atlas
  //   FieldSet. These fields would be 0 in State yy, and they could be part of the norm
  //   computation, so there's no guarantee that xx.norm() == yy.norm() should succeed.
  // - Could compute an Increment dx = xx-yy and check its norm is ~0, but this would make
  //   the State depend on the Increment.
  //
  // So, we skip this test. Instead, go back to a FieldSet and compare FieldSets:
  atlas::FieldSet fset3{};
  yy.toFieldSet(fset3);
  EXPECT(util::compareFieldSets(fset, fset3));
}

// -----------------------------------------------------------------------------
/*! \brief Tests State::geometry() and Geometry copy constructors
 */

template <typename MODEL> void testStateGeometry() {
  typedef StateFixture<MODEL>   Test_;
  typedef oops::Geometry<MODEL> Geometry_;
  typedef oops::State<MODEL>    State_;

  const double norm = Test_::test().normFile;
  const double tol = Test_::test().tolerance;
  const util::DateTime vt(Test_::test().date);

  State_ xx1(Test_::resol(), Test_::test().statefile);

  // get geometry from xx1 and initialize xx2 (xx2 & xx1 should be the same)
  const Geometry_ & geometry = xx1.geometry();
  State_ xx2(geometry, Test_::test().statefile);

  const double norm2 = xx2.norm();
  EXPECT(oops::is_close(norm2, norm, tol));
}

// -----------------------------------------------------------------------------

/*! \brief Interpolation test
 *
 * \details **testStateAnalyticInitialCondition()** tests the creation of an
 * analytic state for a given model.  The conceptual steps are as follows:
 *
 * 1. Initialize the JEDI State object based on idealized analytic formulae
 * 2. Compare norm to precomputed norm
 *
 * Relevant parameters in the **State* section of the config file include
 *
 * * **norm-gen** Normalization test for the generated State
 *
 * \date April, 2018: M. Miesch (JCSDA) adapted a preliminary version in the
 * feature/interp branch
 *
 * \warning Since this model compares the interpolated state values to an exact analytic
 * solution, it requires that the "analytic_init" option be implemented in the model and
 * selected in the "state.state generate" section of the config file.
 */
template <typename MODEL> void testStateAnalyticInitialCondition() {
  typedef StateFixture<MODEL>    Test_;
  typedef oops::State<MODEL>     State_;

  // This creates a State object called xx based on information
  // from the "geometry" and "state test.state generate" sections of
  // the config file and checks its norm

  if (Test_::test().stateGenerate.value() == boost::none ||
      Test_::test().normGeneratedState.value() == boost::none) {
    oops::Log::warning() << "Bypassing Analytical Initial Condition Test";
    return;
  }

  const State_ xx(Test_::resol(), *Test_::test().stateGenerate.value());

  const double norm = *Test_::test().normGeneratedState.value();
  const double tol = Test_::test().tolerance;

  oops::Log::debug() << "xx.norm(): " << std::fixed << std::setprecision(8) << xx.norm()
                     << std::endl;
  oops::Log::debug() << "norm: " << std::fixed << std::setprecision(8) << norm << std::endl;

  EXPECT(oops::is_close(xx.norm(), norm, tol));
}

// -----------------------------------------------------------------------------

/*! \brief Tests of zero and accumul
 *
 * \details testStateZeroAndAccumul tests the folllowing:
 *
 * 1. Call State.zero() and check the resulting norm is indeed zero.
 * 2. Call State.accumul(), passing a multiplication factor and a second State
 *    as arguments, and verify that the resulting norm is as expected.
 *    Several multiplication factors are tested.
 */

template <typename MODEL> void testStateZeroAndAccumul() {
  typedef StateFixture<MODEL>    Test_;
  typedef oops::State<MODEL>     State_;

  State_ xx(Test_::resol(), Test_::test().statefile);
  const double tol = Test_::test().tolerance;

  // Set state xx to zero
  xx.zero();
  EXPECT(xx.norm() == 0.0);

  // Set state xx to various multiples of yy
  const State_ yy(Test_::resol(), Test_::test().statefile);
  const std::vector<double> mults {3.0, 0.0, -3.0};
  for (const auto & mult : mults) {
    xx.zero();
    xx.accumul(mult, yy);
    EXPECT(oops::is_close(xx.norm(), std::abs(mult) * yy.norm(), tol));
  }

  // Ensure that a non-zero state, when acted on with accumul, is not equal to the result
  State_ zz(Test_::resol(), Test_::test().statefile);
  zz.accumul(3.0, yy);
  EXPECT_NOT(oops::is_close(zz.norm(), yy.norm(), tol, 0, oops::TestVerbosity::SILENT));
}

// -----------------------------------------------------------------------------

/*! \brief validTime and updateTime tests
 *
 * \details **testStateDateTime()** tests the validTime and updateTime routines.
 *
 * This is performed by updating the initial state time in two ways:
 * - two lots of one hour,
 * - one lot of two hours.
 *
 * validTime is then used in the comparison of the two times obtained.
 */

template <typename MODEL> void testStateDateTime() {
  typedef StateFixture<MODEL>    Test_;
  typedef oops::State<MODEL>     State_;

  // Configuration to read initial state
  State_ xx(Test_::resol(), Test_::test().statefile);

  // Update the time by two lots of one hour
  const util::Duration onehour(3600);
  xx.updateTime(onehour);
  xx.updateTime(onehour);

  // Create another state
  State_ yy(Test_::resol(), Test_::test().statefile);

  // Update the time of the second state by two hours
  const util::Duration twohours(7200);
  yy.updateTime(twohours);

  EXPECT(xx.validTime() == yy.validTime());

  // Increment the first state's time again and check the times are now not equal
  xx.updateTime(onehour);
  EXPECT_NOT(xx.validTime() == yy.validTime());
}

// -----------------------------------------------------------------------------
/*! \brief Serialize and Deserialize test
 *
 * \details **testStateSerialize()** tests the serialization and deserialization of a state.
*/

template <typename MODEL> void testStateSerialize() {
  typedef StateFixture<MODEL>   Test_;
  typedef oops::State<MODEL>    State_;

  // Configuration to read initial state
  State_ xx(Test_::resol(), Test_::test().statefile);
  const util::Duration tt("PT15H");
  xx.updateTime(tt);

// Create another state
  State_ yy(Test_::resol(), Test_::test().statefile);
  yy.zero();
  yy.accumul(2.5, xx);

  EXPECT(yy.validTime() != xx.validTime());
  if (xx.norm() > 0.0) EXPECT(yy.norm() != xx.norm());

// Test serialize-deserialize
  std::vector<double> vect;
  xx.serialize(vect);
  EXPECT(vect.size() == xx.serialSize());

  size_t index = 0;
  yy.deserialize(vect, index);

  EXPECT(index == xx.serialSize());
  EXPECT(index == yy.serialSize());

  xx.serialize(vect);
  EXPECT(vect.size() == xx.serialSize() * 2);

  if (xx.serialSize() > 0) {  // until all models have implemented serialize
    EXPECT(xx.norm() > 0.0);
    EXPECT(yy.norm() > 0.0);
    EXPECT(yy.validTime() == xx.validTime());
    xx.accumul(-1.0, yy);
    EXPECT(xx.norm() == 0);
  }
}

// -----------------------------------------------------------------------------

/*! \brief Read and write tests
 *
 * \details **testStateReadWrite()** tests reading and writing model state files.
 *
 * The tests are as follows:
 * 1. Read an input file and check it sets a state correctly.
 * 2. Write an output file and read it in again, checking the states are the same.
 */

template <typename MODEL> void testStateReadWrite() {
  typedef StateFixture<MODEL>          Test_;
  typedef oops::State<MODEL>           State_;

  // Configuration to read initial state
  State_ xx(Test_::resol(), Test_::test().statefile);
  const double tol = Test_::test().tolerance;

  // Determine initial state norm
  const double norm = xx.norm();

  // Set state to zero
  xx.zero();

  // Read input file
  xx.read(Test_::test().statefile.value());

  // Check norm has its initial value
  EXPECT(xx.norm() == norm);

  if (Test_::test().writeReadTest.value() != boost::none) {
    const eckit::LocalConfiguration testconf = Test_::test().toConfiguration();
    const eckit::LocalConfiguration rwconf(testconf, "write then read test");
    // Modify state
    const double mult = 2.0;
    xx.accumul(mult, xx);

    // Determine modified state norm
    const double normout = xx.norm();

    // Write modified state to output file
    const eckit::LocalConfiguration wconf(rwconf, "state write");
    xx.write(wconf);

    // Read modified state from output file
    const eckit::LocalConfiguration rconf(rwconf, "state read");
    State_ yy(Test_::resol(), rconf);

    // Check modified state norm has its expected value
    EXPECT(oops::is_close(yy.norm(), normout, tol));

    // Check modified state norm is not equal to the initial norm
    EXPECT_NOT(oops::is_close(norm, normout, tol, 0, oops::TestVerbosity::SILENT));
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class State : public oops::Test {
 public:
  State() {}
  virtual ~State() {StateFixture<MODEL>::reset();}

 private:
  std::string testid() const override {return "test::State<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/State/testStateConstructors")
      { testStateConstructors<MODEL>(); });
    ts.emplace_back(CASE("interface/State/testStateAtlasInterface")
      { testStateAtlasInterface<MODEL>(); });
    ts.emplace_back(CASE("interface/State/testStateGeometry")
      { testStateGeometry<MODEL>(); });
    ts.emplace_back(CASE("interface/State/testStateAnalyticInitialCondition")
      { testStateAnalyticInitialCondition<MODEL>(); });
    ts.emplace_back(CASE("interface/State/testStateZeroAndAccumul")
      { testStateZeroAndAccumul<MODEL>(); });
    ts.emplace_back(CASE("interface/State/testStateDateTime")
      { testStateDateTime<MODEL>(); });
    ts.emplace_back(CASE("interface/State/testStateSerialize")
      { testStateSerialize<MODEL>(); });
    ts.emplace_back(CASE("interface/State/testStateReadWrite")
      { testStateReadWrite<MODEL>(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_STATE_H_
