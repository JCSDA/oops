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

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/PostProcessor.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "test/TestEnvironment.h"

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

    const eckit::LocalConfiguration biasConf(TestEnvironment::config(), "ModelBias");
    bias_.reset(new ModelAux_(*resol_, biasConf));

    const eckit::LocalConfiguration conf(TestEnvironment::config(), "Model");
    model_.reset(new Model_(*resol_, conf));

    const eckit::LocalConfiguration iniConf(TestEnvironment::config(), "State");
    xref_.reset(new State_(*resol_, model_->variables(), iniConf));
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
  EXPECT(Test_::model().timeResolution() > zero);
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

  EXPECT(xx.validTime() == vt);
  EXPECT(xx.norm() == ininorm);

// Recomputing initial norm to make sure nothing bad happened
  EXPECT(Test_::xref().norm() == ininorm);
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

  EXPECT(xx.validTime() == vt);
  EXPECT(oops::is_close(xx.norm(), fnorm, tol));

// Recomputing initial norm to make sure nothing bad happened
  EXPECT(Test_::xref().norm() == ininorm);
}

// =============================================================================

template <typename MODEL>
class Model : public oops::Test {
 public:
  Model() {}
  virtual ~Model() {}
 private:
  std::string testid() const {return "test::Model<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/Model/testModelConstructor")
      { testModelConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/Model/testModelNoForecasat")
      { testModelNoForecast<MODEL>(); });
    ts.emplace_back(CASE("interface/Model/testModelForecast")
      { testModelForecast<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_MODEL_H_
