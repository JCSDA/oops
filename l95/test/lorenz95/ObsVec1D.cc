/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <memory>
#include <string>
#include <vector>


#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/ObsTableView.h"
#include "lorenz95/ObsVec1D.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsVecTestFixture : TestFixture {
 public:
  ObsVecTestFixture() {
    const eckit::LocalConfiguration conf(TestConfig::config(), "observations");
    const util::DateTime bgn(conf.getString("window begin"));
    const util::DateTime end(conf.getString("window end"));
    const eckit::LocalConfiguration otconf(conf, "observation.obs space");
    obstable_.reset(new lorenz95::ObsTableView(otconf, oops::mpi::comm(), bgn, end));
    const std::vector<std::string> vv{"zz"};
  }
  ~ObsVecTestFixture() {}
  std::unique_ptr<lorenz95::ObsTableView> obstable_;
};
// -----------------------------------------------------------------------------
CASE("test_ObsVec1D") {
  ObsVecTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_constructor") {
    std::unique_ptr<lorenz95::ObsVec1D> ov(new lorenz95::ObsVec1D(*fix.obstable_));
    EXPECT(ov.get() != NULL);
    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov)[ii] == 0.0);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_nobs") {
    std::unique_ptr<lorenz95::ObsVec1D> ov(new lorenz95::ObsVec1D(*fix.obstable_));
    EXPECT(ov->nobs() == fix.obstable_->nobs());
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_read") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));
    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov)[ii] != 0.0);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_copy_constructor_copy") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov1(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));
    std::unique_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1));

    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov2)[ii] == (*ov1)[ii]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_classname") {
    std::unique_ptr<lorenz95::ObsVec1D> ov(new lorenz95::ObsVec1D(*fix.obstable_));
    EXPECT(ov->classname() == "lorenz95::ObsVec1D");
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_assignment") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov1(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));
    std::unique_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*fix.obstable_));

    // use the assignment operator to populate the second ObsVec1D object
    *ov2 = *ov1;

    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov2)[ii] == (*ov1)[ii]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_compound_assignment_multiply_double") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov1(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));
    std::unique_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1));

    // create a random double value
    double mult = 7.92;
    *ov2 *= mult;

    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov2)[ii] == (*ov1)[ii] * mult);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_compound_assignment_add") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov1(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));
    std::unique_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1));

    *ov1 += *ov2;

    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov1)[ii] == (*ov2)[ii] + (*ov2)[ii]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_compound_assignment_subtract") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov1(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));
    std::unique_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1));

    *ov1 -= *ov2;

    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov1)[ii] == 0.0);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_compound_assignment_multiply") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov1(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));
    std::unique_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1));

    *ov1 *= *ov2;

    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov1)[ii] == (*ov2)[ii] * (*ov2)[ii]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_compound_assignment_divide") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov1(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));
    std::unique_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1));

    *ov1 /= *ov2;

    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov1)[ii] == 1.0);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_zero") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));

    ov->zero();

    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov)[ii] == 0.0);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_axpy") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov1(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));
    std::unique_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1));

    double mult = 2.00;

    ov1->axpy(mult, *ov2);

    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov1)[ii] == (*ov2)[ii] + mult * (*ov2)[ii]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_invert") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov1(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));
    std::unique_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1));

    ov1->invert();

    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov1)[ii] == 1.0 / (*ov2)[ii]);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_random") {
    std::unique_ptr<lorenz95::ObsVec1D> ov(new lorenz95::ObsVec1D(*fix.obstable_));

    ov->random();

    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      EXPECT((*ov)[ii] != 0);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_ObsVec1D_dot_product_with") {
    std::unique_ptr<lorenz95::ObsVec1D>
      ov1(new lorenz95::ObsVec1D(*fix.obstable_, "ObsValue"));
    std::unique_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1));

    double result = ov1->dot_product_with(*ov2);

    double check = 0.0;
    for (unsigned int ii = 0; ii < fix.obstable_->nobs(); ++ii) {
      check += (*ov2)[ii] * (*ov2)[ii];
    }

    EXPECT(oops::is_close(result, check, 1.0e-10));
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
