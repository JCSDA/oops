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
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/IgnoreOtherParameters.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

/// Options used by testStateReadWrite().
template <typename MODEL>
class StateWriteReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StateWriteReadParameters, Parameters)

 public:
  typedef oops::State<MODEL>                State_;
  typedef typename State_::Parameters_      StateParameters_;
  typedef typename State_::WriteParameters_ StateWriteParameters_;

  /// Options used by the code writing the state to a file.
  oops::RequiredParameter<StateWriteParameters_> write{"state write", this};
  /// Options used by the code reading the state back in.
  oops::RequiredParameter<StateParameters_> read{"state read", this};
};

// -----------------------------------------------------------------------------

/// Configuration of the state test.
template <typename MODEL>
class StateTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StateTestParameters, Parameters)

 public:
  typedef oops::State<MODEL>                 State_;
  typedef StateWriteReadParameters<MODEL>    StateWriteReadParameters_;
  typedef typename State_::Parameters_       StateParameters_;

  /// Relative tolerance of norm comparisons.
  oops::RequiredParameter<double> tolerance{"tolerance", this};
  /// Validity time for states loaded from a file and generated on the fly.
  oops::RequiredParameter<util::DateTime> date{"date", this};

  /// Configuration of the state loaded from a file.
  oops::RequiredParameter<StateParameters_> statefile{"statefile", this};
  /// Expected norm of the state loaded from a file.
  oops::RequiredParameter<double> normFile{"norm file", this};

  /// Configuration of the state generated on the fly.
  oops::OptionalParameter<eckit::LocalConfiguration> stateGenerate{"state generate", this};
  /// Expected norm of the state generated on the fly.
  ///
  /// This option must be present if `state generate` is.
  oops::OptionalParameter<double> normGeneratedState{"norm generated state", this};

  oops::OptionalParameter<StateWriteReadParameters_> writeReadTest{"write then read test", this};
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
  xx.read(Test_::test().statefile);

  // Check norm has its initial value
  EXPECT(xx.norm() == norm);

  if (Test_::test().writeReadTest.value() != boost::none) {
    // Modify state
    const double mult = 2.0;
    xx.accumul(mult, xx);

    // Determine modified state norm
    const double normout = xx.norm();

    // Write modified state to output file
    xx.write(Test_::test().writeReadTest.value()->write);

    // Read modified state from output file
    State_ yy(Test_::resol(), Test_::test().writeReadTest.value()->read);

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
