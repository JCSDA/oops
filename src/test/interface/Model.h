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
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Geometry.h"
#include "oops/base/Model.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "test/TestEnvironment.h"

namespace test {

// =================================================================================================

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
  static void reset() {
    getInstance().xref_.reset();
    getInstance().bias_.reset();
    getInstance().model_.reset();
    getInstance().resol_.reset();
    getInstance().test_.reset();
  }

 private:
  static ModelFixture<MODEL>& getInstance() {
    static ModelFixture<MODEL> theModelFixture;
    return theModelFixture;
  }

  ModelFixture<MODEL>() {
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "model test"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::world()));

    const eckit::LocalConfiguration biasConf(TestEnvironment::config(), "model aux control");
    bias_.reset(new ModelAux_(*resol_, biasConf));

    const eckit::LocalConfiguration conf(TestEnvironment::config(), "model");
    model_.reset(new Model_(*resol_, conf));

    const eckit::LocalConfiguration iniConf(TestEnvironment::config(), "initial condition");
    xref_.reset(new State_(*resol_, iniConf));
  }

  ~ModelFixture<MODEL>() {}

  std::unique_ptr<const eckit::LocalConfiguration>  test_;
  std::unique_ptr<const Geometry_>     resol_;
  std::unique_ptr<const State_>        xref_;
  std::unique_ptr<const ModelAux_>     bias_;
  std::unique_ptr<const Model_>        model_;
};

// =================================================================================================

/// \brief tests constructor, timeResolution() method and print method
template <typename MODEL> void testModelConstructor() {
  typedef ModelFixture<MODEL>   Test_;
  const util::Duration zero(0);
  EXPECT(Test_::model().timeResolution() > zero);
  oops::Log::test() << "Testing Model: " << Test_::model() << std::endl;
}

// -------------------------------------------------------------------------------------------------

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

// -------------------------------------------------------------------------------------------------

template <typename MODEL> void testModelForecast() {
  typedef ModelFixture<MODEL>   Test_;
  typedef oops::State<MODEL>    State_;

  const double fnorm = Test_::test().getDouble("final norm");
  const double tol = Test_::test().getDouble("tolerance");
  const util::Duration len(Test_::test().getString("forecast length"));

  const double ininorm = Test_::xref().norm();
  State_ xx(Test_::xref());
  const util::DateTime vt(xx.validTime()+len);

  oops::PostProcessor<State_> post;

  Test_::model().forecast(xx, Test_::bias(), len, post);

  EXPECT(xx.validTime() == vt);

  oops::Log::debug() << "xx.norm(): " << std::fixed << std::setprecision(8) << xx.norm()
                     << std::endl;
  oops::Log::debug() << "fnorm: " << std::fixed << std::setprecision(8) << fnorm << std::endl;

  EXPECT(oops::is_close(xx.norm(), fnorm, tol));

// Recomputing initial norm to make sure nothing bad happened
  EXPECT(Test_::xref().norm() == ininorm);
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL> void testModelReForecast() {
  typedef ModelFixture<MODEL>   Test_;
  typedef oops::State<MODEL>    State_;

  const bool testreforecast = Test_::test().getBool("test reforecast", true);

  if (testreforecast) {
    // Forecast duration
    const util::Duration len(Test_::test().getString("forecast length"));

    const double ininorm = Test_::xref().norm();

    // Initial states
    State_ xx1(Test_::xref());
    State_ xx2(Test_::xref());

    // Time at forecast end
    const util::DateTime vt(xx1.validTime()+len);

    oops::PostProcessor<State_> post;

    // Forecast 1
    Test_::model().forecast(xx1, Test_::bias(), len, post);

    // Forecast 2
    Test_::model().forecast(xx2, Test_::bias(), len, post);

    // Check forecasts ran to expected time
    EXPECT(xx1.validTime() == vt);
    EXPECT(xx2.validTime() == vt);

    // Print the final norms
    oops::Log::debug() << "xx1.norm(): " << std::fixed << std::setprecision(8) << xx1.norm()
                       << std::endl;
    oops::Log::debug() << "xx2.norm(): " << std::fixed << std::setprecision(8) << xx2.norm()
                       << std::endl;

    // Pass or fail condition
    EXPECT(xx1.norm() == xx2.norm());

    // Recomputing initial norm to make sure nothing bad happened
    EXPECT(Test_::xref().norm() == ininorm);
  } else {
    // Dummy test
    EXPECT(0.0 == 0.0);
  }
}

// =================================================================================================

template <typename MODEL>
class Model : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~Model() {ModelFixture<MODEL>::reset();}
 private:
  std::string testid() const override {return "test::Model<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/Model/testModelConstructor")
      { testModelConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/Model/testModelNoForecast")
      { testModelNoForecast<MODEL>(); });
    ts.emplace_back(CASE("interface/Model/testModelForecast")
      { testModelForecast<MODEL>(); });
    ts.emplace_back(CASE("interface/Model/testModelReForecast")
      { testModelReForecast<MODEL>(); });
  }

  void clear() const override {}
};

// =================================================================================================

}  // namespace test

#endif  // TEST_INTERFACE_MODEL_H_
