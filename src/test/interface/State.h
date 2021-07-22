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
#include "oops/base/Variables.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

template <typename MODEL> class StateFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>       Geometry_;
  typedef oops::State<MODEL>          State_;

 public:
  static const eckit::Configuration & test()  {return *getInstance().test_;}
  static const Geometry_    & resol() {return *getInstance().resol_;}
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
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "state test"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::world()));
  }

  ~StateFixture<MODEL>() {}

  std::unique_ptr<const eckit::LocalConfiguration>  test_;
  std::unique_ptr<Geometry_>     resol_;
};

// -----------------------------------------------------------------------------
/// \brief tests constructors and print method
template <typename MODEL> void testStateConstructors() {
  typedef StateFixture<MODEL>   Test_;
  typedef oops::State<MODEL>    State_;

  const double norm = Test_::test().getDouble("norm file");
  const double tol = Test_::test().getDouble("tolerance");
  const util::DateTime vt(Test_::test().getString("date"));

// Test main constructor
  const eckit::LocalConfiguration conf(Test_::test(), "statefile");
  std::unique_ptr<State_> xx1(new State_(Test_::resol(), conf));

  EXPECT(xx1.get());
  oops::Log::test() << "Printing State from yaml: " << *xx1 << std::endl;
  const double norm1 = xx1->norm();
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

  const double norm = Test_::test().getDouble("norm file");
  const double tol = Test_::test().getDouble("tolerance");
  const util::DateTime vt(Test_::test().getString("date"));

  const eckit::LocalConfiguration conf(Test_::test(), "statefile");
  State_ xx1(Test_::resol(), conf);

  // get geometry from xx1 and initialize xx2 (xx2 & xx1 should be the same)
  const Geometry_ & geometry = xx1.geometry();
  State_ xx2(geometry, conf);

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

  if (!Test_::test().has("state generate")) {
    oops::Log::warning() << "Bypassing Analytical Initial Condition Test";
    return;
  }

  const eckit::LocalConfiguration confgen(Test_::test(), "state generate");
  const State_ xx(Test_::resol(), confgen);

  const double norm = Test_::test().getDouble("norm generated state");
  const double tol = Test_::test().getDouble("tolerance");

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

  const eckit::LocalConfiguration conf(Test_::test(), "statefile");
  State_ xx(Test_::resol(), conf);
  const double tol = Test_::test().getDouble("tolerance");

  // Set state xx to zero
  xx.zero();
  EXPECT(xx.norm() == 0.0);

  // Set state xx to various multiples of yy
  const State_ yy(Test_::resol(), conf);
  const std::vector<double> mults {3.0, 0.0, -3.0};
  for (const auto & mult : mults) {
    xx.zero();
    xx.accumul(mult, yy);
    EXPECT(oops::is_close(xx.norm(), std::abs(mult) * yy.norm(), tol));
  }

  // Ensure that a non-zero state, when acted on with accumul, is not equal to the result
  State_ zz(Test_::resol(), conf);
  zz.accumul(3.0, yy);
  EXPECT_NOT(oops::is_close(zz.norm(), yy.norm(), tol, 0, oops::TestVerbosity::SILENT));
}

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
  const eckit::LocalConfiguration conf(Test_::test(), "statefile");
  State_ xx(Test_::resol(), conf);

  // Update the time by two lots of one hour
  const util::Duration onehour(3600);
  xx.updateTime(onehour);
  xx.updateTime(onehour);

  // Create another state
  State_ yy(Test_::resol(), conf);

  // Update the time of the second state by two hours
  const util::Duration twohours(7200);
  yy.updateTime(twohours);

  EXPECT(xx.validTime() == yy.validTime());

  // Increment the first state's time again and check the times are now not equal
  xx.updateTime(onehour);
  EXPECT_NOT(xx.validTime() == yy.validTime());
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
  typedef StateFixture<MODEL>    Test_;
  typedef oops::State<MODEL>     State_;

  // Configuration to read initial state
  const eckit::LocalConfiguration conf(Test_::test(), "statefile");
  State_ xx(Test_::resol(), conf);
  const double tol = Test_::test().getDouble("tolerance");

  // Determine initial state norm
  const double norm = xx.norm();

  // Set state to zero
  xx.zero();

  // Read input file
  xx.read(conf);

  // Check norm has its initial value
  EXPECT(xx.norm() == norm);

  if (Test_::test().has("statefileout")) {
    // Modify state
    const double mult = 2.0;
    xx.accumul(mult, xx);

    // Determine modified state norm
    const double normout = xx.norm();

    // Configuration to read and write output state
    const eckit::LocalConfiguration confout(Test_::test(), "statefileout");

    // Write modified state to output file
    xx.write(confout);

    // Read modifed state from output file
    State_ yy(Test_::resol(), confout);

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
    ts.emplace_back(CASE("interface/State/testStateReadWrite")
      { testStateReadWrite<MODEL>(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_STATE_H_
