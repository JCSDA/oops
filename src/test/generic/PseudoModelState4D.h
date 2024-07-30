/*
 * (C) Copyright 2021- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_GENERIC_PSEUDOMODELSTATE4D_H_
#define TEST_GENERIC_PSEUDOMODELSTATE4D_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Geometry.h"
#include "oops/base/Model.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State4D.h"
#include "oops/generic/PseudoModelState4D.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Duration.h"
#include "test/TestEnvironment.h"

namespace test {

/// \brief Tests that PseudoModelState4D, run to 0th, 1st, 2nd, etc
/// state produces the same results as the Model specified in yaml.
template <typename MODEL> void testPseudoModelState4D() {
  typedef oops::Geometry<MODEL>        Geometry_;
  typedef oops::Model<MODEL>           Model_;
  typedef oops::ModelAuxControl<MODEL> ModelAux_;
  typedef oops::State<MODEL>           State_;
  typedef oops::State4D<MODEL>         State4D_;
  typedef oops::ModelBase<MODEL>       ModelBase_;
  typedef oops::PseudoModelState4D<MODEL> PseudoModelState4D_;

  // Setup geometry, model bias, and initial conditions
  const eckit::LocalConfiguration geometryconf(TestEnvironment::config(), "geometry");
  Geometry_ geometry(geometryconf, oops::mpi::world());
  const eckit::LocalConfiguration biasconf(TestEnvironment::config(), "model aux control");
  ModelAux_ modelaux(geometry, biasconf);
  const eckit::LocalConfiguration initconf(TestEnvironment::config(), "initial condition");
  State_ xinit(geometry, initconf);

  // set up Model specified in yaml
  const eckit::LocalConfiguration modelconf(TestEnvironment::config(), "model");
  Model_ model1(geometry, modelconf);
  oops::Log::test() << "Model: " << model1 << std::endl;

  // set up PseudoModelState4D with State4D from yaml
  const eckit::LocalConfiguration state4dconf(TestEnvironment::config(),
                                              "pseudo model with 4D state");
  State4D_ state4d(geometry, state4dconf);
  std::unique_ptr<ModelBase_> pseudomodel(new PseudoModelState4D_(state4d));
  Model_ model2(std::move(pseudomodel));
  oops::Log::test() << "PseudoModelState4D: " << model2 << std::endl;

  // check that the two models have the same time resolutions
  const util::Duration step = model1.timeResolution();
  EXPECT(step == model2.timeResolution());
  const size_t nsteps = state4d.size();

  // run forecasts with 0, step, 2*step, ... lengths and compare the results
  for (size_t jstep = 0; jstep < nsteps; ++jstep) {
    State_ x1(xinit);
    State_ x2(xinit);
    oops::PostProcessor<State_> post;
    oops::Log::test() << "Running Model for " << jstep << " steps" << std::endl;
    model1.forecast(x1, modelaux, step * jstep, post);
    oops::Log::test() << "Running PseudoModelState4D for " << jstep << " steps" << std::endl;
    model2.forecast(x2, modelaux, step * jstep, post);
    EXPECT(x1.norm() == x2.norm());
    EXPECT(x1.norm() == state4d[jstep].norm());
  }
}

template <typename MODEL>
class PseudoModelState4D : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~PseudoModelState4D() = default;
 private:
  std::string testid() const override {return "test::PseudoModelState4D<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("generic/PseudoModelState4D/testPseudoModelState4D")
      { testPseudoModelState4D<MODEL>(); });
  }

  void clear() const override {}
};

}  // namespace test

#endif  // TEST_GENERIC_PSEUDOMODELSTATE4D_H_
