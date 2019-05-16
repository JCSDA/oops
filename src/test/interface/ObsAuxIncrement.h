/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_OBSAUXINCREMENT_H_
#define TEST_INTERFACE_OBSAUXINCREMENT_H_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxCovariance.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "test/TestEnvironment.h"

namespace test {

// =============================================================================

template <typename MODEL> class ObsAuxIncrementFixture : private boost::noncopyable {
  typedef oops::ObsAuxCovariance<MODEL> Covariance_;
  typedef oops::ObsAuxControl<MODEL>    ObsAux_;
  typedef oops::ObsAuxIncrement<MODEL>  AuxIncr_;

 public:
  static const eckit::Configuration & config()     {return *getInstance().conf_;}
  static const Covariance_  & covariance() {return *getInstance().covar_;}

 private:
  static ObsAuxIncrementFixture<MODEL>& getInstance() {
    static ObsAuxIncrementFixture<MODEL> theObsAuxIncrementFixture;
    return theObsAuxIncrementFixture;
  }

  ObsAuxIncrementFixture<MODEL>() {
    std::vector<eckit::LocalConfiguration> osconf;
    TestEnvironment::config().get("Observations.ObsTypes", osconf);
    conf_.reset(new eckit::LocalConfiguration(osconf[0]));
    covar_.reset(new Covariance_(*conf_));
  }

  ~ObsAuxIncrementFixture<MODEL>() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
  boost::scoped_ptr<const Covariance_>  covar_;
};

// =============================================================================

template <typename MODEL> void testObsAuxIncrementConstructor() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx(Test_::config());

  EXPECT(dx.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementCopyConstructor() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

  AuxIncr_ dx2(dx1);
  EXPECT(dx2.norm() > 0.0);
  EXPECT(dx2.norm() == dx1.norm());

// Check that the copy is equal to the original
  dx2 -= dx1;
  EXPECT(dx2.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementChangeRes() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

  AuxIncr_ dx2(dx1, Test_::config());
  EXPECT(dx2.norm() > 0.0);
  EXPECT(dx2.norm() == dx1.norm());

// Check that the copy is equal to the original
  dx2 -= dx1;
  EXPECT(dx2.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementTriangle() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx2);

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

template <typename MODEL> void testObsAuxIncrementOpPlusEq() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(dx1);

// test *= and +=
  dx2 += dx1;
  dx1 *= 2.0;

  dx2 -= dx1;
  EXPECT(dx2.norm() < 1e-8);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementDotProduct() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);
  AuxIncr_ dx2(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx2);

// test symmetry of dot product
  double zz1 = dot_product(dx1, dx2);
  double zz2 = dot_product(dx2, dx1);

  EXPECT(zz1 == zz2);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementZero() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx);
  EXPECT(dx.norm() > 0.0);

// test zero
  dx->zero();
  EXPECT(dx.norm() == 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL> void testObsAuxIncrementAxpy() {
  typedef ObsAuxIncrementFixture<MODEL>   Test_;
  typedef oops::ObsAuxIncrement<MODEL>    AuxIncr_;

  AuxIncr_ dx1(Test_::config());
  ObsAuxIncrementFixture<MODEL>::covariance().randomize(dx1);

// test axpy
  AuxIncr_ dx2(dx1);
  dx2.axpy(2.0, dx1);

  dx2 -= dx1;
  dx2 -= dx1;
  dx2 -= dx1;

  EXPECT(dx2.norm() < 1e-8);
}

// =============================================================================

template <typename MODEL>
class ObsAuxIncrement : public oops::Test {
 public:
  ObsAuxIncrement() {}
  virtual ~ObsAuxIncrement() {}

 private:
  std::string testid() const {return "test::ObsAuxIncrement<" + MODEL::name() + ">";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementConstructor")
      { testObsAuxIncrementConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementCopyConstructor")
      { testObsAuxIncrementCopyConstructor<MODEL>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementChangeRes")
      { testObsAuxIncrementChangeRes<MODEL>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementTriangle")
      { testObsAuxIncrementTriangle<MODEL>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementOpPlusEq")
      { testObsAuxIncrementOpPlusEq<MODEL>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementDotProduct")
      { testObsAuxIncrementDotProduct<MODEL>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementAxpy")
      { testObsAuxIncrementAxpy<MODEL>(); });
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSAUXINCREMENT_H_
