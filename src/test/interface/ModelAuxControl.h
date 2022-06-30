/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_MODELAUXCONTROL_H_
#define TEST_INTERFACE_MODELAUXCONTROL_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/parameters/IgnoreOtherParameters.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

/// \brief Parameters loaded from the input YAML file and used by this test.
template <typename MODEL>
class TestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TestParameters, Parameters)

 public:
  typedef typename oops::Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef typename oops::ModelAuxControl<MODEL>::Parameters_ ModelAuxControlParameters_;

  /// \brief Group of parameters controlling the tested model's geometry.
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};
  /// \brief Group of parameters controlling the tested implementation of the ModelAuxControl
  /// interface.
  oops::RequiredParameter<ModelAuxControlParameters_> modelAuxControl{"model aux control", this};
  /// \brief Don't treat the presence of other parameter groups as an error (this makes it
  /// possible to reuse a single YAML file in tests of implementations of multiple oops interfaces).
  oops::IgnoreOtherParameters ignoreOthers{this};
};

// -----------------------------------------------------------------------------
template <typename MODEL> class ModelAuxControlFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>        Geometry_;
  typedef oops::ModelAuxControl<MODEL> ModelAux_;
  typedef TestParameters<MODEL>        TestParameters_;

 public:
  static const typename ModelAux_::Parameters_ & parameters() {return getInstance().parameters_;}
  static const Geometry_    & resol()  {return *getInstance().resol_;}

 private:
  static ModelAuxControlFixture<MODEL>& getInstance() {
    static ModelAuxControlFixture<MODEL> theModelAuxControlFixture;
    return theModelAuxControlFixture;
  }

  ModelAuxControlFixture() {
    TestParameters_ parameters;
    parameters.validateAndDeserialize(TestEnvironment::config());

    resol_.reset(new Geometry_(parameters.geometry, oops::mpi::world()));
    parameters_ = parameters.modelAuxControl;
  }

  ~ModelAuxControlFixture() {}

  typename ModelAux_::Parameters_ parameters_;
  std::unique_ptr<Geometry_>     resol_;
};

// -----------------------------------------------------------------------------
/// \brief testing constructor and print method
template <typename MODEL> void testConstructor() {
  typedef ModelAuxControlFixture<MODEL>   Test_;
  typedef oops::ModelAuxControl<MODEL>    ModelAux_;

  std::unique_ptr<ModelAux_> bias(new ModelAux_(Test_::resol(), Test_::parameters()));
  EXPECT(bias.get());
  oops::Log::test() << "Testing ModelAuxControl: " << *bias << std::endl;
  bias.reset();
  EXPECT(!bias.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef ModelAuxControlFixture<MODEL>   Test_;
  typedef oops::ModelAuxControl<MODEL>    ModelAux_;

  std::unique_ptr<ModelAux_> bias(new ModelAux_(Test_::resol(), Test_::parameters()));

  std::unique_ptr<ModelAux_> other(new ModelAux_(*bias));
  EXPECT(other.get());

  other.reset();
  EXPECT(!other.get());

  EXPECT(bias.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testChangeRes() {
  typedef ModelAuxControlFixture<MODEL>   Test_;
  typedef oops::ModelAuxControl<MODEL>    ModelAux_;

  std::unique_ptr<ModelAux_> bias(new ModelAux_(Test_::resol(), Test_::parameters()));

  std::unique_ptr<ModelAux_> other(new ModelAux_(Test_::resol(), *bias));
  EXPECT(other.get());

  other.reset();
  EXPECT(!other.get());

  EXPECT(bias.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class ModelAuxControl : public oops::Test {
 public:
  ModelAuxControl() {}
  virtual ~ModelAuxControl() {}
 private:
  std::string testid() const override {return "test::ModelAuxControl<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ModelAuxControl/testConstructor")
      { testConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ModelAuxControl/testCopyConstructor")
      { testCopyConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ModelAuxControl/testChangeRes")
      { testChangeRes<MODEL>(); });
  }

  void clear() const override {}
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_MODELAUXCONTROL_H_
