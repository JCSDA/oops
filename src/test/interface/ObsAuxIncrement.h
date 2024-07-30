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
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxCovariance.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// =============================================================================

template <typename OBS> class ObsAuxIncrementFixture : private boost::noncopyable {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsAuxCovariance<OBS> Covariance_;

 public:
  static const Covariance_    & covariance(const size_t ii) {return *getInstance().covar_.at(ii);}

 private:
  static ObsAuxIncrementFixture<OBS>& getInstance() {
    static ObsAuxIncrementFixture<OBS> theObsAuxIncrementFixture;
    return theObsAuxIncrementFixture;
  }

  ObsAuxIncrementFixture<OBS>() {
    for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
      eckit::LocalConfiguration biasconf = Test_::config(jj).getSubConfiguration("obs bias");
      std::shared_ptr<Covariance_> tmp(new Covariance_(Test_::obspace()[jj], biasconf));
      covar_.push_back(tmp);
    }
  }

  ~ObsAuxIncrementFixture<OBS>() {}

  std::vector<std::shared_ptr<Covariance_> > covar_;
};

// =============================================================================
/// \brief test constructor and print method
template <typename OBS> void testObsAuxIncrementConstructor() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsAuxIncrement<OBS>    AuxIncr_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration biasconf = Test_::config(jj).getSubConfiguration("obs bias");
    AuxIncr_ dx(Test_::obspace()[jj], biasconf);
    oops::Log::info() << "Printing zero ObsAuxIncrement: " << dx << std::endl;
    EXPECT(dx.norm() == 0.0);
  }
}

// -----------------------------------------------------------------------------
/// Tests copy-constructor (with option of allocating, but not copying the data)
template <typename OBS> void testObsAuxIncrementCopyConstructor() {
  typedef ObsTestsFixture<OBS>        Test_;
  typedef oops::ObsAuxIncrement<OBS>  AuxIncr_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    if (Test_::config(jj).has("obs bias")) {
      eckit::LocalConfiguration biasconf = Test_::config(jj).getSubConfiguration("obs bias");
      AuxIncr_ dx1(Test_::obspace()[jj], biasconf);
      ObsAuxIncrementFixture<OBS>::covariance(jj).randomize(dx1);
      oops::Log::info() << "Printing random ObsAuxIncrement: " << dx1 << std::endl;
      /// Test that creating new increment without copying data works
      AuxIncr_ dxempty(dx1, false);
      EXPECT_EQUAL(dxempty.norm(), 0.0);
      /// Test that creating new increment with copying data works
      AuxIncr_ dx2(dx1);
      EXPECT(dx2.norm() > 0.0);
      EXPECT(dx2.norm() == dx1.norm());

      /// Test that the copy is equal to the original
      dx2 -= dx1;
      EXPECT(dx2.norm() == 0.0);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testObsAuxIncrementTriangle() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef ObsAuxIncrementFixture<OBS>   AuxTest_;
  typedef oops::ObsAuxIncrement<OBS>    AuxIncr_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    if (Test_::config(jj).has("obs bias")) {
      eckit::LocalConfiguration biasconf = Test_::config(jj).getSubConfiguration("obs bias");
      AuxIncr_ dx1(Test_::obspace()[jj], biasconf);
      AuxTest_::covariance(jj).randomize(dx1);
      AuxIncr_ dx2(Test_::obspace()[jj], biasconf);
      ObsAuxIncrementFixture<OBS>::covariance(jj).randomize(dx2);

//    test triangle inequality
      double dot1 = dx1.norm();
      EXPECT(dot1 > 0.0);

      double dot2 = dx2.norm();
      EXPECT(dot2 > 0.0);

      dx2 += dx1;
      double dot3 = dx2.norm();
      EXPECT(dot3 > 0.0);

      EXPECT(dot3 <= dot1 + dot2);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testObsAuxIncrementOpPlusEq() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef ObsAuxIncrementFixture<OBS>   AuxTest_;
  typedef oops::ObsAuxIncrement<OBS>    AuxIncr_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration biasconf = Test_::config(jj).getSubConfiguration("obs bias");
    AuxIncr_ dx1(Test_::obspace()[jj], biasconf);
    AuxTest_::covariance(jj).randomize(dx1);
    AuxIncr_ dx2(dx1);

//  test *= and +=
    dx2 += dx1;
    dx1 *= 2.0;

    dx2 -= dx1;
    EXPECT(dx2.norm() < 1e-8);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testObsAuxIncrementDotProduct() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef ObsAuxIncrementFixture<OBS>   AuxTest_;
  typedef oops::ObsAuxIncrement<OBS>    AuxIncr_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration biasconf = Test_::config(jj).getSubConfiguration("obs bias");
    AuxIncr_ dx1(Test_::obspace()[jj], biasconf);
    AuxTest_::covariance(jj).randomize(dx1);
    AuxIncr_ dx2(Test_::obspace()[jj], biasconf);
    AuxTest_::covariance(jj).randomize(dx2);

//  test symmetry of dot product
    double zz1 = dot_product(dx1, dx2);
    double zz2 = dot_product(dx2, dx1);

    EXPECT(zz1 == zz2);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testObsAuxIncrementZero() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef ObsAuxIncrementFixture<OBS>   AuxTest_;
  typedef oops::ObsAuxIncrement<OBS>    AuxIncr_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration biasconf = Test_::config(jj).getSubConfiguration("obs bias");
    AuxIncr_ dx(Test_::obspace()[jj], biasconf);
    AuxTest_::covariance(jj).randomize(dx);
    EXPECT(dx.norm() > 0.0);

//  test zero
    dx->zero();
    EXPECT(dx.norm() == 0.0);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testObsAuxIncrementAxpy() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef ObsAuxIncrementFixture<OBS>   AuxTest_;
  typedef oops::ObsAuxIncrement<OBS>    AuxIncr_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration biasconf = Test_::config(jj).getSubConfiguration("obs bias");
    AuxIncr_ dx1(Test_::obspace()[jj], biasconf);
    AuxTest_::covariance(jj).randomize(dx1);

//  test axpy
    AuxIncr_ dx2(dx1);
    dx2.axpy(2.0, dx1);

    dx2 -= dx1;
    dx2 -= dx1;
    dx2 -= dx1;

    EXPECT(dx2.norm() < 1e-8);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testObsAuxIncrementSerializeDeserialize() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef ObsAuxIncrementFixture<OBS>   AuxTest_;
  typedef oops::ObsAuxIncrement<OBS>    AuxIncr_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    eckit::LocalConfiguration biasconf = Test_::config(jj).getSubConfiguration("obs bias");
    AuxIncr_ dx1(Test_::obspace()[jj], biasconf);
    AuxIncr_ dx2(Test_::obspace()[jj], biasconf);
    dx2.zero();

    AuxTest_::covariance(jj).randomize(dx1);

    EXPECT(dx1.norm() != dx2.norm());

    // test serialize / deserialize
    std::vector<double> vect;
    dx1.serialize(vect);
    EXPECT(vect.size() == dx1.serialSize());

    size_t index = 0;
    dx2.deserialize(vect, index);
    EXPECT(index == dx1.serialSize());
    EXPECT(index == dx2.serialSize());

    dx1.serialize(vect);
    EXPECT(vect.size() == dx1.serialSize() * 2);
    EXPECT(dx1.norm() > 0.0);
    EXPECT(dx2.norm() > 0.0);

    dx2 -= dx1;
    EXPECT(dx2.norm() == 0.0);
  }
}

// =============================================================================

template <typename OBS>
class ObsAuxIncrement : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;

 public:
  using oops::Test::Test;
  virtual ~ObsAuxIncrement() {}

 private:
  std::string testid() const override {return "test::ObsAuxIncrement<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementConstructor")
      { testObsAuxIncrementConstructor<OBS>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementCopyConstructor")
      { testObsAuxIncrementCopyConstructor<OBS>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementTriangle")
      { testObsAuxIncrementTriangle<OBS>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementOpPlusEq")
      { testObsAuxIncrementOpPlusEq<OBS>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementDotProduct")
      { testObsAuxIncrementDotProduct<OBS>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementAxpy")
      { testObsAuxIncrementAxpy<OBS>(); });
    ts.emplace_back(CASE("interface/ObsAuxIncrement/testObsAuxIncrementSerializeDeserialize")
      { testObsAuxIncrementAxpy<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSAUXINCREMENT_H_
