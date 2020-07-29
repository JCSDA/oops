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
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
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

 private:
  static StateFixture<MODEL>& getInstance() {
    static StateFixture<MODEL> theStateFixture;
    return theStateFixture;
  }

  StateFixture<MODEL>() {
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "state test"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::comm()));
  }

  ~StateFixture<MODEL>() {}

  std::unique_ptr<const eckit::LocalConfiguration>  test_;
  std::unique_ptr<Geometry_>     resol_;
};

// -----------------------------------------------------------------------------

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

template <typename MODEL>
class State : public oops::Test {
 public:
  State() {}
  virtual ~State() {}
 private:
  std::string testid() const {return "test::State<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/State/testStateConstructors")
      { testStateConstructors<MODEL>(); });
    ts.emplace_back(CASE("interface/State/testStateAnalyticInitialCondition")
      { testStateAnalyticInitialCondition<MODEL>(); });
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_STATE_H_
