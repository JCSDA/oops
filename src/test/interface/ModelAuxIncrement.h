/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_MODELAUXINCREMENT_H_
#define TEST_INTERFACE_MODELAUXINCREMENT_H_

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
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxCovariance.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "test/TestEnvironment.h"

namespace test {

// =============================================================================

template <typename MODEL> class ModelAuxIncrementFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>           Geometry_;
  typedef oops::ModelAuxCovariance<MODEL> Covariance_;
  typedef oops::ModelAuxControl<MODEL>    ModelAux_;
  typedef oops::ModelAuxIncrement<MODEL>  AuxIncr_;

 public:
  static const eckit::Configuration & config()     {return *getInstance().conf_;}
  static const Covariance_  & covariance() {return *getInstance().covar_;}
  static const Geometry_    & resol()      {return *getInstance().resol_;}

 private:
  static ModelAuxIncrementFixture<MODEL>& getInstance() {
    static ModelAuxIncrementFixture<MODEL> theModelAuxIncrementFixture;
    return theModelAuxIncrementFixture;
  }

  ModelAuxIncrementFixture<MODEL>() {
//  Setup a geometry
    const eckit::LocalConfiguration resolConfig(TestEnvironment::config(), "geometry");
    resol_.reset(new Geometry_(resolConfig, oops::mpi::world()));

//  Setup a covariance matrix
    conf_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "model aux error"));
    covar_.reset(new Covariance_(*conf_, *resol_));
  }

  ~ModelAuxIncrementFixture<MODEL>() {}

  std::unique_ptr<const eckit::LocalConfiguration> conf_;
  std::unique_ptr<const Geometry_>    resol_;
  std::unique_ptr<const Covariance_>  covar_;
};

// =============================================================================
/// \brief tests constructor and print method
template <typename MODEL> void testModelAuxIncrementConstructor() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx(Test_::resol(), Test_::config());
  oops::Log::test() << "Testing ModelAuxIncrement: " << dx << std::endl;
  EXPECT(dx.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementCopyConstructor() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

  AuxIncr_ dx2(dx1);
  EXPECT(dx2.norm() > 0.0);
  EXPECT(dx2.norm() == dx1.norm());

// Check that the copy is equal to the original
  dx2 -= dx1;
  EXPECT(dx2.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementChangeRes() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

  AuxIncr_ dx2(dx1, Test_::config());
  EXPECT(dx2.norm() > 0.0);
  EXPECT(dx2.norm() == dx1.norm());

// Check that the copy is equal to the original
  dx2 -= dx1;
  EXPECT(dx2.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementTriangle() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx2);

// test triangle inequality
  double dot1 = dx1.norm();
  EXPECT(dot1 > 0.0);

  double dot2 = dx2.norm();
  EXPECT(dot2 > 0.0);

  dx2 += dx1;
  double dot3 = dx2.norm();
  EXPECT(dot3 > 0.0);

  EXPECT(dot3 <= dot1 + dot2);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementOpPlusEq() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(dx1);

// test *= and +=
  dx2 += dx1;
  dx1 *= 2.0;

  dx2 -= dx1;
  EXPECT(dx2.norm() < 1e-8);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementDotProduct() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx2);

// test symmetry of dot product
  double zz1 = dot_product(dx1, dx2);
  double zz2 = dot_product(dx2, dx1);

  EXPECT(zz1 == zz2);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementZero() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx);
  EXPECT(dx.norm() > 0.0);

// test zero
  dx->zero();
  EXPECT(dx.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testModelAuxIncrementAxpy() {
  typedef ModelAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ModelAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::resol(), Test_::config());
  ModelAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

// test axpy
  AuxIncr_ dx2(dx1);
  dx2.axpy(2.0, dx1);

  dx2 -= dx1;
  dx2 -= dx1;
  dx2 -= dx1;

  EXPECT(dx2.norm() < 1e-8);
}

// =============================================================================

template <typename MODEL> class ModelAuxIncrement : public oops::Test {
 public:
  ModelAuxIncrement() {}
  virtual ~ModelAuxIncrement() {}

 private:
  std::string testid() const override {return "test::ModelAuxIncrement<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();
    ts.emplace_back(CASE("interface/ModelAuxIncrement/testModelAuxIncrementConstructor")
      { testModelAuxIncrementConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ModelAuxIncrement/testModelAuxIncrementCopyConstructor")
      { testModelAuxIncrementCopyConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ModelAuxIncrement/testModelAuxIncrementChangeRes")
      { testModelAuxIncrementChangeRes<MODEL>(); });
    ts.emplace_back(CASE("interface/ModelAuxIncrement/testModelAuxIncrementTriangle")
      { testModelAuxIncrementTriangle<MODEL>(); });
    ts.emplace_back(CASE("interface/ModelAuxIncrement/testModelAuxIncrementOpPlusEq")
      { testModelAuxIncrementOpPlusEq<MODEL>(); });
    ts.emplace_back(CASE("interface/ModelAuxIncrement/testModelAuxIncrementDotProduct")
      { testModelAuxIncrementDotProduct<MODEL>(); });
    ts.emplace_back(CASE("interface/ModelAuxIncrement/testModelAuxIncrementAxpy")
      { testModelAuxIncrementAxpy<MODEL>(); });
  }

  void clear() const override {}
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_MODELAUXINCREMENT_H_
