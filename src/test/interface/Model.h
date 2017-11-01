/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_MODEL_H_
#define TEST_INTERFACE_MODEL_H_

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
#include "oops/base/PostProcessor.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "test/TestEnvironment.h"
#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace test {

// =============================================================================

template <typename MODEL> class ModelFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>        Geometry_;
  typedef oops::Model<MODEL>           Model_;
  typedef oops::ModelAuxControl<MODEL> ModelAux_;
  typedef oops::State<MODEL>           State_;

 public:
  static const eckit::Configuration & test()  {return *getInstance().test_;}
  static const Geometry_    & resol() {return *getInstance().resol_;}
  static const State_       & xref()  {return *getInstance().xref_;}
  static const ModelAux_    & bias()  {return *getInstance().bias_;}
  static const Model_       & model() {return *getInstance().model_;}

 private:
  static ModelFixture<MODEL>& getInstance() {
    static ModelFixture<MODEL> theModelFixture;
    return theModelFixture;
  }

  ModelFixture<MODEL>() {
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "ModelTest"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));

    const eckit::LocalConfiguration iniConf(TestEnvironment::config(), "State");
    xref_.reset(new State_(*resol_, iniConf));

    const eckit::LocalConfiguration biasConf(TestEnvironment::config(), "ModelBias");
    bias_.reset(new ModelAux_(*resol_, biasConf));

    const eckit::LocalConfiguration conf(TestEnvironment::config(), "Model");
    model_.reset(new Model_(*resol_, conf));
  }

  ~ModelFixture<MODEL>() {}

  boost::scoped_ptr<const eckit::LocalConfiguration>  test_;
  boost::scoped_ptr<const Geometry_>     resol_;
  boost::scoped_ptr<const State_>        xref_;
  boost::scoped_ptr<const ModelAux_>     bias_;
  boost::scoped_ptr<const Model_>        model_;
};

// =============================================================================

template <typename MODEL> void testModelConstructor() {
  typedef ModelFixture<MODEL>   Test_;

  const util::Duration zero(0);
  BOOST_CHECK(Test_::model().timeResolution() > zero);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelNoForecast() {
  typedef ModelFixture<MODEL>   Test_;
  typedef oops::State<MODEL>    State_;

  const double ininorm = Test_::xref().norm();
  State_ xx(Test_::xref());
  const util::DateTime vt(xx.validTime());

  const util::Duration zero(0);
  oops::PostProcessor<State_> post;

  Test_::model().forecast(xx, Test_::bias(), zero, post);

  BOOST_CHECK_EQUAL(xx.validTime(), vt);
  BOOST_CHECK_EQUAL(xx.norm(), ininorm);

// Recomputing initial norm to make sure nothing bad happened
  BOOST_CHECK_EQUAL(Test_::xref().norm(), ininorm);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelForecast() {
  typedef ModelFixture<MODEL>   Test_;
  typedef oops::State<MODEL>    State_;

  const double fnorm = Test_::test().getDouble("finalnorm");
  const double tol = Test_::test().getDouble("tolerance");
  const util::Duration len(Test_::test().getString("fclength"));

  const double ininorm = Test_::xref().norm();
  State_ xx(Test_::xref());
  const util::DateTime vt(xx.validTime()+len);

  oops::PostProcessor<State_> post;

  Test_::model().forecast(xx, Test_::bias(), len, post);

  BOOST_CHECK_EQUAL(xx.validTime(), vt);
  BOOST_CHECK_CLOSE(xx.norm(), fnorm, tol);

// Recomputing initial norm to make sure nothing bad happened
  BOOST_CHECK_EQUAL(Test_::xref().norm(), ininorm);
}

// =============================================================================

template <typename MODEL> class Model : public oops::Test {
 public:
  Model() {}
  virtual ~Model() {}
 private:
  std::string testid() const {return "test::Model<" + MODEL::name() + ">";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("interface/Model");

    ts->add(BOOST_TEST_CASE(&testModelConstructor<MODEL>));
    ts->add(BOOST_TEST_CASE(&testModelNoForecast<MODEL>));
    ts->add(BOOST_TEST_CASE(&testModelForecast<MODEL>));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_MODEL_H_
