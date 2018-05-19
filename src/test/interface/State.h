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
#include <string>
#include <vector>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
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
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "StateTest"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));
  }

  ~StateFixture<MODEL>() {}

  boost::scoped_ptr<const eckit::LocalConfiguration>  test_;
  boost::scoped_ptr<Geometry_>     resol_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testStateConstructors() {
  typedef StateFixture<MODEL>   Test_;
  typedef oops::State<MODEL>    State_;

  const double norm = Test_::test().getDouble("norm-file");
  const double tol = Test_::test().getDouble("tolerance");
  const util::DateTime vt(Test_::test().getString("date"));

// Test main constructor
  const eckit::LocalConfiguration conf(Test_::test(), "StateFile");
  boost::scoped_ptr<State_> xx1(new State_(Test_::resol(), conf));

  BOOST_CHECK(xx1.get());
  const double norm1 = xx1->norm();
  BOOST_CHECK_CLOSE(norm1, norm, tol);
  BOOST_CHECK_EQUAL(xx1->validTime(), vt);

// Test copy constructor
  boost::scoped_ptr<State_> xx2(new State_(*xx1));
  BOOST_CHECK(xx2.get());
  BOOST_CHECK_CLOSE(xx2->norm(), norm, tol);
  BOOST_CHECK_EQUAL(xx2->validTime(), vt);

// Destruct copy
  xx2.reset();
  BOOST_CHECK(!xx2.get());

// Recompute initial norm to make sure nothing bad happened
  const double norm2 = xx1->norm();
  BOOST_CHECK_EQUAL(norm1, norm2);
}

// -----------------------------------------------------------------------------
/*! \brief Interpolation test
 *
 * \details **testStateInterpolation()** tests the interpolation for a given
 * model.  The conceptual steps are as follows:
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
 * user in the "StateTest" section of the config file in addition to a
 * randomly-generated list of **Nrandom** random locations.  Nrandom is also
 * specified by the user in the "StateTest" section of the config file, as is the
 * (nondimensional) tolerence level (**interp_tolerance**) to be used for the tests.
 *
 * Relevant parameters in the **State* section of the config file include
 *
 * * **norm-gen** Normalization test for the generated State
 * * **interp_tolerance** tolerance for the interpolation test
 *
 * \date April, 2018: M. Miesch (JCSDA) adapted a preliminary version in the
 * feature/interp branch
 *
 * \warning Since this model compares the interpolated state values to an exact analytic
 * solution, it requires that the "analytic_init" option be implemented in the model and
 * selected in the "State.StateGenerate" section of the config file.
 */

template <typename MODEL> void testStateInterpolation() {
  typedef StateFixture<MODEL>    Test_;
  typedef oops::State<MODEL>     State_;
  typedef oops::Locations<MODEL> Locations_;
  typedef oops::GeoVaLs<MODEL>   GeoVaLs_;

  // This creates a State object called xx based on information
  // from the "Geometry" and "StateTest.StateGenerate" sections of
  // the config file and checks its norm

  const eckit::LocalConfiguration confgen(Test_::test(), "StateGenerate");
  const State_ xx(Test_::resol(), confgen);
  const double norm = Test_::test().getDouble("norm-gen");
  const double tol = Test_::test().getDouble("tolerance");
  BOOST_CHECK_CLOSE(xx.norm(), norm, tol);

  // Now extract the user-defined locations from the "StateTest.Locations"
  // section of the config file and use it to define a Locations object
  // The user can optionally also request Nrandom random locations
  const eckit::LocalConfiguration confloc(Test_::test(), "Locations");
  const Locations_ locs(confloc);

  // Extract the user-defined list of variables to interpolate,
  // also from the "StateGenerate" section of the config file, and
  // use this to define a Variables object
  const oops::Variables vars(confgen);

  // Now create a GeoVaLs object from locs and vars
  GeoVaLs_ gval(locs, vars);

  // ...and execute the interpolation
  xx.getValues(locs, vars, gval);

  // Now create another GeoVaLs object that contains the exact
  // analytic solutions
  GeoVaLs_ ref(locs, vars, confgen);

  // Compute the difference between the interpolated and exact values
  gval -= ref;

  // Compute the normalized error
  gval.abs();
  gval /= ref;

  // And check to see if the errors are within specified tolerance
  double interp_tol = Test_::test().getDouble("interp_tolerance");
  BOOST_CHECK_SMALL(gval.norm(), interp_tol);

  // Each MODEL should define an appropriate GeoVaLs print() method that
  // writes information about the GeoVaLs object to the oops::Log::debug()
  // output stream in order to help with debugging in the event that the
  // interpolation test does not pass.  It is helpful if this print()
  // method at least prints out the location and variable where the
  // error is largest.

  oops::Log::debug() << "TestStateInterpolation() Normalized Error: "
                     << std::endl << gval << std::endl;

  // Also write the locations used for the test to the debug output stream
  oops::Log::debug() << "TestStateInterpolation() Locations: "
                     << std::endl << locs << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL> class State : public oops::Test {
 public:
  State() {}
  virtual ~State() {}
 private:
  std::string testid() const {return "test::State<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/State");

    ts->add(BOOST_TEST_CASE(&testStateConstructors<MODEL>));
    ts->add(BOOST_TEST_CASE(&testStateInterpolation<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_STATE_H_
