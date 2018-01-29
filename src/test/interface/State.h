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

#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"
#include "util/Logger.h"
#include "util/dot_product.h"

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

template <typename MODEL> void testStateInterpolation() {
  typedef StateFixture<MODEL>     Test_;
  typedef oops::State<MODEL>      State_;
  typedef oops::Locations<MODEL>  Locations_;
  typedef oops::GeoVaLs<MODEL>    GeoVaLs_;

  const eckit::LocalConfiguration confs(Test_::test(), "StateGenerate");
  const State_ xx(Test_::resol(), confs);
  const double norm = Test_::test().getDouble("norm-gen");
  const double tol = Test_::test().getDouble("tolerance");
  BOOST_CHECK_CLOSE(xx.norm(), norm, tol);

  const eckit::LocalConfiguration confl(Test_::test(), "Locations");
  const Locations_ locs(confl);

  const eckit::LocalConfiguration confv(Test_::test(), "Variables");
  const oops::Variables vars(confv);

  GeoVaLs_ gval(locs, vars);

  xx.interpolate(locs, vars, gval);

  std::vector<double> values;
  Test_::test().get("values", values);
  if (values.size() > 0) {
    const double zz = std::sqrt(dot_product(gval, gval));
    const double ref = values[0];
    BOOST_CHECK_CLOSE(zz, ref, 0.5);
  }
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
