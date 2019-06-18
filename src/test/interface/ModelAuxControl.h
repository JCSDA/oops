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
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/runs/Test.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
template <typename MODEL> class ModelAuxControlFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>        Geometry_;
  typedef oops::ModelAuxControl<MODEL> ModelAux_;

 public:
  static const eckit::Configuration & config() {return *getInstance().conf_;}
  static const Geometry_    & resol()  {return *getInstance().resol_;}

 private:
  static ModelAuxControlFixture<MODEL>& getInstance() {
    static ModelAuxControlFixture<MODEL> theModelAuxControlFixture;
    return theModelAuxControlFixture;
  }

  ModelAuxControlFixture() {
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "ModelBias"));

    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "Geometry");
    resol_.reset(new Geometry_(resolConfig));
  }

  ~ModelAuxControlFixture() {}

  std::unique_ptr<const eckit::LocalConfiguration>  conf_;
  std::unique_ptr<Geometry_>     resol_;
};

// -----------------------------------------------------------------------------

template <typename MODEL> void testConstructor() {
  typedef ModelAuxControlFixture<MODEL>   Test_;
  typedef oops::ModelAuxControl<MODEL>    ModelAux_;

  std::unique_ptr<ModelAux_> bias(new ModelAux_(Test_::resol(), Test_::config()));
  EXPECT(bias.get());

  bias.reset();
  EXPECT(!bias.get());
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testCopyConstructor() {
  typedef ModelAuxControlFixture<MODEL>   Test_;
  typedef oops::ModelAuxControl<MODEL>    ModelAux_;

  std::unique_ptr<ModelAux_> bias(new ModelAux_(Test_::resol(), Test_::config()));

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

  std::unique_ptr<ModelAux_> bias(new ModelAux_(Test_::resol(), Test_::config()));

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
  std::string testid() const {return "test::ModelAuxControl<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ModelAuxControl/testConstructor")
      { testConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ModelAuxControl/testCopyConstructor")
      { testCopyConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ModelAuxControl/testChangeRes")
      { testChangeRes<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_MODELAUXCONTROL_H_
