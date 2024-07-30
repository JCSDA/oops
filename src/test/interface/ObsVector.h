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
#include <utility>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/testing/Test.h"
#include "oops/base/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsDataVector.h"
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
    oops::Log::info() << "Printing zero ObsVector: " << *ov << std::endl;
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
    oops::Log::info() << "Printing random ObsVector: " << *ov << std::endl;

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

/// Test the constructor taking a std::unique_ptr<OBS::ObsVector>.
template <typename OBS> void testWrappingConstructor() {
  typedef ObsTestsFixture<OBS>  Test_;
  typedef oops::ObsVector<OBS>  ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    ObsVector_ ov(Test_::obspace()[jj]);

    ov.random();
    oops::Log::info() << "Printing random ObsVector: " << ov << std::endl;

    ObsVector_ other(std::make_unique<typename OBS::ObsVector_>(ov->obsvector()),
                     Test_::obspace()[jj]);

    const double ov2 = dot_product(ov, ov);
    const double other2 = dot_product(other, other);

    EXPECT(ov2 == other2);
  }
}

// -----------------------------------------------------------------------------
/// Test that:
/// - dot_product of a random vector with self is non-zero
/// - dot_product of a random vector with a zero vector is zero
/// - dot_product of a zero vector with self is zero
/// - dot_product of a vector of ones with self is equal to nobs
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

    ObsVector_ ov3(ov1);
    ov3.ones();

    EXPECT(dot_product(ov3, ov3) == ov3.nobs());
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
/// \brief Tests ObsVector::mask, ObsVector::packEigen and
///        ObsVector::packEigenSize methods.
/// \details Tests that:
/// - mask of "nothing to mask" applied to ObsVector doesn't change its size
///   and content;
/// - mask of either "mask all" (if "mask variable" isn't specified in yaml), or
///   from the file applied to ObsVector changes its size.
/// - linear algebra operations with ObsVector that were masked out produce
///   ObsVectors that have the same number of obs masked out.
/// - size returned by packEigenSize is consistent with size of Eigen Vector
///   returned by packEigen, and is the same as reference value for each MPI
///   task.
template <typename OBS> void testMask() {
  typedef ObsTestsFixture<OBS>           Test_;
  typedef oops::ObsDataVector<OBS, int>  ObsDataVectorInt_;
  typedef oops::ObsDataVector<OBS, float>  ObsDataVector_;
  typedef oops::ObsSpace<OBS>            ObsSpace_;
  typedef oops::ObsVector<OBS>           ObsVector_;

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
    const ObsSpace_ & obspace = Test_::obspace()[jj];

    ObsVector_ reference(obspace);
    reference.random();
    oops::Log::info() << "ObsVector before masking: " << reference << std::endl;

    const size_t nobs_all = Test_::config(jj).getInt("reference global nobs");
    EXPECT_EQUAL(reference.nobs(), nobs_all);
    EXPECT(nobs_all > 0);

    /// apply empty mask, check that vector is the same
    ObsVector_ unsetmask(obspace);
    unsetmask.zero();
    ObsVector_ with_unsetmask(reference);
    with_unsetmask.mask(unsetmask);
    oops::Log::info() << "ObsVector masked with all-zero-mask: " << with_unsetmask << std::endl;
    EXPECT_EQUAL(with_unsetmask.nobs(), nobs_all);
    with_unsetmask -= reference;
    EXPECT_EQUAL(with_unsetmask.rms(), 0.0);

    /// apply non-empty mask, check number of observations
    std::string maskvarname;
    /// by default (else statement below), apply mask that masks out everything
    size_t nobs_after_mask = 0;
    std::vector<size_t> nobs_after_mask_local(Test_::comm().size(), 0);
    /// if mask variable is available, apply mask from file and read reference number of masked obs
    if (Test_::config(jj).has("mask variable")) {
      maskvarname = Test_::config(jj).getString("mask variable");
      nobs_after_mask = Test_::config(jj).getUnsigned("reference global masked nobs");
      nobs_after_mask_local = Test_::config(jj).getUnsignedVector("reference local masked nobs");
      // check that specified mask masks out something
      EXPECT_NOT_EQUAL(nobs_after_mask, nobs_all);
      // check that "reference local masked nobs" are defined for all MPI tasks
      EXPECT_EQUAL(Test_::comm().size(), nobs_after_mask_local.size());
    /// if mask variable is unavailable, apply mask with all ones
    } else {
      // Hack for mask with ones: use ObsVector set to ones, write to ObsSpace,
      // then read as ObsDataVector.
      maskvarname = "set_mask";
      ObsVector_ tmp(obspace);
      tmp.ones();
      tmp.save(maskvarname);
    }
    // emulate the flow in the applications: ObsDataVector_ obs errors
    // get masked with ObsDataVectorInt_ QC flags; then copied into ObsVector
    // obs errors, and used as a mask for another ObsVector (e.g. H(x)).
    ObsDataVectorInt_ qcmask(obspace, obspace.assimvariables(), maskvarname);
    ObsDataVector_ masked_withqc(obspace, obspace.assimvariables());
    ObsVector_ mask(obspace);
    masked_withqc.mask(qcmask);
    mask = masked_withqc;

    ObsVector_ with_mask(reference);
    with_mask.mask(mask);
    oops::Log::info() << "ObsVector masked with " << maskvarname << " mask: " <<
                         with_mask << std::endl;
    EXPECT_EQUAL(with_mask.nobs(), nobs_after_mask);

    /// Test that various linear algebra operations with masked out ObsVector
    /// produce masked out ObsVector
    /// Test invert()
    ObsVector_ test(with_mask);
    test.invert();
    EXPECT_EQUAL(test.nobs(), nobs_after_mask);

    /// Test *=(float)
    test *= 2.0;
    EXPECT_EQUAL(test.nobs(), nobs_after_mask);

    /// Test +=(ObsVector)
    test.random();
    EXPECT_EQUAL(test.nobs(), nobs_all);
    test += with_mask;
    EXPECT_EQUAL(test.nobs(), nobs_after_mask);
    test.random();
    EXPECT_EQUAL(test.nobs(), nobs_all);
    with_mask += test;
    EXPECT_EQUAL(with_mask.nobs(), nobs_after_mask);

    /// Test -=(ObsVector)
    EXPECT_EQUAL(test.nobs(), nobs_all);
    test -= with_mask;
    EXPECT_EQUAL(test.nobs(), nobs_after_mask);
    test.random();
    EXPECT_EQUAL(test.nobs(), nobs_all);
    with_mask -= test;
    EXPECT_EQUAL(with_mask.nobs(), nobs_after_mask);

    /// Test *=(ObsVector)
    EXPECT_EQUAL(test.nobs(), nobs_all);
    test *= with_mask;
    EXPECT_EQUAL(test.nobs(), nobs_after_mask);
    test.random();
    EXPECT_EQUAL(test.nobs(), nobs_all);
    with_mask *= test;
    EXPECT_EQUAL(with_mask.nobs(), nobs_after_mask);

    /// Test /=(ObsVector)
    EXPECT_EQUAL(test.nobs(), nobs_all);
    test /= with_mask;
    EXPECT_EQUAL(test.nobs(), nobs_after_mask);
    test.random();
    EXPECT_EQUAL(test.nobs(), nobs_all);
    with_mask /= test;
    EXPECT_EQUAL(with_mask.nobs(), nobs_after_mask);

    /// Test axpy
    EXPECT_EQUAL(test.nobs(), nobs_all);
    test.axpy(2.0, with_mask);
    EXPECT_EQUAL(test.nobs(), nobs_after_mask);
    test.random();
    EXPECT_EQUAL(test.nobs(), nobs_all);
    with_mask.axpy(2.0, test);
    EXPECT_EQUAL(with_mask.nobs(), nobs_after_mask);

    /// test packEigen
    // create maskvec ObsVector - mask with missing values
    ObsVector_ maskvec = test;
    maskvec.ones();
    maskvec.mask(mask);
    // randomize the vector, and call packEigen with maskvec
    test.random();
    Eigen::VectorXd with_mask_vec = test.packEigen(maskvec);
    // check that the size of returned Eigen Vector is consistent with size
    // returned by packEigenSize()
    EXPECT_EQUAL(static_cast<size_t>(with_mask_vec.size()), test.packEigenSize(maskvec));
    oops::Log::debug() << "Local number of masked observations is: " <<
                          with_mask_vec.size() << std::endl;
    // check that the size is consistent with reference for this MPI task
    EXPECT_EQUAL(static_cast<size_t>(with_mask_vec.size()),
                 nobs_after_mask_local[Test_::comm().rank()]);
  }
}
// -----------------------------------------------------------------------------

template <typename OBS>
class ObsVector : public oops::Test {
  typedef ObsTestsFixture<OBS> Test_;

 public:
  using oops::Test::Test;
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
    ts.emplace_back(CASE("interface/ObsVector/testMask")
      { testMask<OBS>(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// =============================================================================

}  // namespace test

#endif  // TEST_INTERFACE_OBSVECTOR_H_
