/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_MODELAUXCOVARIANCE_H_
#define TEST_INTERFACE_MODELAUXCOVARIANCE_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Geometry.h"
#include "oops/interface/ModelAuxCovariance.h"
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
  typedef typename oops::ModelAuxCovariance<MODEL>::Parameters_ ModelAuxCovarianceParameters_;

  /// \brief Group of parameters controlling the tested model's geometry.
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};
  /// \brief Group of parameters controlling the tested implementation of the ModelAuxCovariance
  /// interface.
  oops::RequiredParameter<ModelAuxCovarianceParameters_> modelAuxError{"model aux error", this};
  /// \brief Don't treat the presence of other parameter groups as an error (this makes it
  /// possible to reuse a single YAML file in tests of implementations of multiple oops interfaces).
  oops::IgnoreOtherParameters ignoreOthers{this};
};

// -----------------------------------------------------------------------------
template <typename MODEL> class ModelAuxCovarianceFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>           Geometry_;
  typedef oops::ModelAuxCovariance<MODEL> Covariance_;
  typedef TestParameters<MODEL>           TestParameters_;

 public:
  static const typename Covariance_::Parameters_ & parameters() {return getInstance().parameters_;}
  static const Geometry_    & resol()  {return *getInstance().resol_;}

 private:
  static ModelAuxCovarianceFixture<MODEL>& getInstance() {
    static ModelAuxCovarianceFixture<MODEL> theModelAuxCovarianceFixture;
    return theModelAuxCovarianceFixture;
  }

  ModelAuxCovarianceFixture() {
    TestParameters_ parameters;
    parameters.validateAndDeserialize(TestEnvironment::config());

    parameters_ = parameters.modelAuxError;

    resol_.reset(new Geometry_(parameters.geometry, oops::mpi::world()));
  }

  ~ModelAuxCovarianceFixture() {}

  typename Covariance_::Parameters_ parameters_;
  std::unique_ptr<Geometry_>     resol_;
};

// -----------------------------------------------------------------------------
/// \brief tests constructor and print method
template <typename MODEL> void testConstructor() {
  typedef ModelAuxCovarianceFixture<MODEL>   Test_;
  typedef oops::ModelAuxCovariance<MODEL>    Covariance_;

  std::unique_ptr<Covariance_> cov(new Covariance_(Test_::parameters(), Test_::resol()));
  EXPECT(cov.get());
  oops::Log::test() << "Testing ModelAuxCovariance: " << *cov << std::endl;
  cov.reset();
  EXPECT(!cov.get());
}

// -----------------------------------------------------------------------------
//  void linearize(const ModelAuxControl_ &, const Geometry_ &);
//  void multiply(const ModelAuxIncrement_ &, ModelAuxIncrement_ &) const;
//  void invMult(const ModelAuxIncrement_ &, ModelAuxIncrement_ &) const;
//  void randomize(ModelAuxIncrement_ &) const;
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

template <typename MODEL>
class ModelAuxCovariance : public oops::Test {
 public:
  ModelAuxCovariance() {}
  virtual ~ModelAuxCovariance() {}
 private:
  std::string testid() const override {return "test::ModelAuxCovariance<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ModelAuxCovariance/testConstructor")
      { testConstructor<MODEL>(); });
  }

  void clear() const override {}
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_MODELAUXCOVARIANCE_H_
