/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_OBSVECTOR_H_
#define TEST_INTERFACE_OBSVECTOR_H_

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/testing/Test.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsVector.h"
#include "oops/runs/Test.h"
#include "oops/util/dot_product.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// \brief tests constructor and print method
template <typename OBS> void testConstructor() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsVector<OBS>  ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    std::unique_ptr<ObsVector_> ov(new ObsVector_(Test_::obspace()[jj]));
    EXPECT(ov.get());
    oops::Log::test() << "Printing zero ObsVector: " << *ov << std::endl;
    ov.reset();
    EXPECT(!ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testCopyConstructor() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsVector<OBS>  ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    std::unique_ptr<ObsVector_> ov(new ObsVector_(Test_::obspace()[jj]));

    ov->random();
    oops::Log::test() << "Printing random ObsVector: " << *ov << std::endl;

    std::unique_ptr<ObsVector_> other(new ObsVector_(*ov));
    EXPECT(other.get());

    const double ov2 = dot_product(*ov, *ov);
    const double other2 = dot_product(*other, *other);

    EXPECT(ov2 == other2);

    other.reset();
    EXPECT(!other.get());

    EXPECT(ov.get());
  }
}

// -----------------------------------------------------------------------------

template <typename OBS> void testNotZero() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsVector<OBS>  ObsVector_;
  const double zero = 0.0;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ ov1(Test_::obspace()[jj]);
    ov1.random();

    const double zz = dot_product(ov1, ov1);
    EXPECT(zz > zero);

    ObsVector_ ov2(ov1);
    ov2.zero();

    EXPECT(dot_product(ov2, ov1) == zero);
    EXPECT(dot_product(ov2, ov2) == zero);
  }
}
// -----------------------------------------------------------------------------

template <typename OBS> void testLinearAlgebra() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsVector<OBS>  ObsVector_;
  const double tolerance = 1.0e-8;
  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ ov1(Test_::obspace()[jj]);
    ov1.random();

    // test *=, += and -=
    ObsVector_ ov2(ov1);
    ov2 += ov1;
    ov1 *= 2.0;
    ov2 -= ov1;
    EXPECT(dot_product(ov2, ov2) < tolerance);

    // test =
    ObsVector_ ov3(ov1);
    ov2 = ov1;
    ov2 -= ov3;
    EXPECT(dot_product(ov2, ov2) < tolerance);

    // test *=(const ObsVector &) and /=(const ObsVector &)
    ov2 = ov1;
    ov2 *= ov1;
    ov2 /= ov1;
    ov2 -= ov1;
    EXPECT(dot_product(ov2, ov2) < tolerance);

    // test axpy
    ov2 = ov1;
    ov3 = ov1;
    ov2.axpy(2.0, ov1);
    ov3 *= 3;
    ov2 -= ov3;
    EXPECT(dot_product(ov2, ov2) < tolerance);

    // test invert()
    ov2 = ov1;
    ov2.invert();
    ov2 *= ov1;
    EXPECT(std::abs(dot_product(ov2, ov2) - ov2.nobs()) < tolerance);
  }
}

// -----------------------------------------------------------------------------
template <typename OBS> void testReadWrite() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsVector<OBS>  ObsVector_;
  const double tolerance = 1.0e-8;
  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ ov1(Test_::obspace()[jj]);
    ov1.random();

    ov1.save("test");
    ObsVector_ ov2(Test_::obspace()[jj], "test");
    ov2 -= ov1;
    EXPECT(dot_product(ov2, ov2) < tolerance);
  }
}
// -----------------------------------------------------------------------------
template <typename OBS> void testPackEigen() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsVector<OBS>  ObsVector_;
  const double tolerance = 1.0e-8;
  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ ov1(Test_::obspace()[jj]);
    ov1.random();
    double rms1 = ov1.rms();

    Eigen::VectorXd vec = ov1.packEigen();
    EXPECT(vec.size() == ov1.nobs());

    double rms2 = sqrt(vec.squaredNorm() / ov1.nobs());
    EXPECT(std::abs(rms1-rms2) < tolerance);
  }
}
// -----------------------------------------------------------------------------

template <typename OBS>
class ObsVector : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;

 public:
  ObsVector() {}
  virtual ~ObsVector() {}

 private:
  std::string testid() const override {return "test::ObsVector<" + OBS::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/ObsVector/testConstructor")
      { testConstructor<OBS>(); });
    ts.emplace_back(CASE("interface/ObsVector/testCopyConstructor")
      { testCopyConstructor<OBS>(); });
    ts.emplace_back(CASE("interface/ObsVector/testNotZero")
      { testNotZero<OBS>(); });
    ts.emplace_back(CASE("interface/ObsVector/testLinearAlgebra")
      { testLinearAlgebra<OBS>(); });
    ts.emplace_back(CASE("interface/ObsVector/testReadWrite")
      { testReadWrite<OBS>(); });
    ts.emplace_back(CASE("interface/ObsVector/testPackEigen")
      { testPackEigen<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSVECTOR_H_
