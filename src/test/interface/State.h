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

#include <iostream>
#include <string>
#include <cmath>

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/runs/Test.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"

namespace test {

// =============================================================================

template <typename MODEL> class StateFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>       Geometry_;
  typedef oops::State<MODEL>          State_;
  typedef oops::Variables<MODEL>      Variables_;

 public:
  static const eckit::Configuration & test()  {return *getInstance().test_;}
  static const Geometry_    & resol() {return *getInstance().resol_;}
  static const Variables_   & vars()  {return *getInstance().vars_;}
  static const State_       & xref()  {return *getInstance().xref_;}
  static const double       & norm()  {return getInstance().refnorm_;}

 private:
  static StateFixture<MODEL>& getInstance() {
    static StateFixture<MODEL> theStateFixture;
    return theStateFixture;
  }

  StateFixture<MODEL>() {
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "StateTest"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));

    const eckit::LocalConfiguration varConfig(TestEnvironment::config(), "Variables");
    vars_.reset(new Variables_(varConfig));

    const eckit::LocalConfiguration conf(TestEnvironment::config(), "State");
    xref_.reset(new State_(*resol_, conf));
    refnorm_ = xref_->norm();
  }

  ~StateFixture<MODEL>() {}

  boost::scoped_ptr<const eckit::LocalConfiguration>  test_;
  boost::scoped_ptr<Geometry_>     resol_;
  boost::scoped_ptr<Variables_>    vars_;
  boost::scoped_ptr<State_>        xref_;
  double                           refnorm_;
};

// =============================================================================

// template <typename MODEL> void testStateConstructor() {
//   typedef StateFixture<MODEL>   Test_;
//   typedef oops::State<MODEL>    State_;
//   const util::DateTime tt(Test_::test().getString("date"));
//   State_ xx(Test_::resol(), Test_::vars(), tt);
// }

// -----------------------------------------------------------------------------

template <typename MODEL> void testStateConstructor() {
  typedef StateFixture<MODEL>   Test_;

  BOOST_CHECK_EQUAL(Test_::xref().norm(), Test_::norm());

  double norm = Test_::test().getDouble("norm");
  double tol = Test_::test().getDouble("tolerance");
  BOOST_CHECK_CLOSE(Test_::xref().norm(), norm, tol);

  const util::DateTime vt(Test_::test().getString("date"));
  BOOST_CHECK_EQUAL(Test_::xref().validTime(), vt);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testStateCopyConstructor() {
  typedef StateFixture<MODEL>   Test_;
  typedef oops::State<MODEL>    State_;

  boost::scoped_ptr<State_> xx(new State_(Test_::xref()));
  BOOST_CHECK(xx.get());
  BOOST_CHECK_EQUAL(xx->norm(), Test_::norm());

  const util::DateTime vt(Test_::test().getString("date"));
  BOOST_CHECK_EQUAL(xx->validTime(), vt);

  xx.reset();
  BOOST_CHECK(!xx.get());

// Recomputing initial norm to make sure nothing bad happened
  BOOST_CHECK_EQUAL(Test_::xref().norm(), Test_::norm());
}

// =============================================================================

template <typename MODEL> class State : public oops::Test {
 public:
  State() {}
  virtual ~State() {}
 private:
  std::string testid() const {return "test::State<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/State");

    ts->add(BOOST_TEST_CASE(&testStateConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testStateCopyConstructor<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_STATE_H_
