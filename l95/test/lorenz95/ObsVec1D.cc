/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "./TestConfig.h"
#include "lorenz95/ObsVec1D.h"
#include "lorenz95/ObsTable.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class ObsVecTestFixture : TestFixture {
 public:
  ObsVecTestFixture() {
    const eckit::LocalConfiguration conf(TestConfig::config(), "Observations");
    const util::DateTime bgn(conf.getString("window_begin"));
    const util::DateTime end(conf.getString("window_end"));
    const eckit::LocalConfiguration otconf(conf, "Observation");
    obstable_.reset(new lorenz95::ObsTable(otconf, bgn, end));
  }
  ~ObsVecTestFixture() {}
  boost::scoped_ptr<lorenz95::ObsTable> obstable_;
};
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
BOOST_FIXTURE_TEST_SUITE(test_ObsVec1D, ObsVecTestFixture)
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_constructor) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov(new lorenz95::ObsVec1D(*obstable_));
    BOOST_CHECK(ov.get() != NULL);
    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov)(ii), 0.0);
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_nobs) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov(new lorenz95::ObsVec1D(*obstable_));
    BOOST_CHECK_EQUAL(ov->size(), obstable_->nobs());
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_read) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov(new lorenz95::ObsVec1D(*obstable_));
    ov->read("ObsVal");
    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK((*ov)(ii) != 0.0);
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_copy_constructor_copy) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov1(new lorenz95::ObsVec1D(*obstable_));
    ov1->read("ObsVal");
    boost::scoped_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1, true));

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov2)(ii), (*ov1)(ii));
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_copy_constructor_no_copy) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov1(new lorenz95::ObsVec1D(*obstable_));
    ov1->read("ObsVal");
    boost::scoped_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1, false));

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov2)(ii), 0.0);
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_classname) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov(new lorenz95::ObsVec1D(*obstable_));
    BOOST_CHECK_EQUAL(ov->classname(), "lorenz95::ObsVec1D");
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_assignment) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov1(new lorenz95::ObsVec1D(*obstable_));
    ov1->read("ObsVal");
    boost::scoped_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*obstable_));

    // use the assignment operator to populate the second ObsVec1D object
    *ov2 = *ov1;

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov2)(ii), (*ov1)(ii));
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_compound_assignment_multiply_double) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov1(new lorenz95::ObsVec1D(*obstable_));
    ov1->read("ObsVal");
    boost::scoped_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1, true));

    // create a random double value
    double mult = 7.92;
    *ov2 *= mult;

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov2)(ii), (*ov1)(ii) * mult);
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_compound_assignment_add) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov1(new lorenz95::ObsVec1D(*obstable_));
    ov1->read("ObsVal");
    boost::scoped_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1, true));

    *ov1 += *ov2;

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov1)(ii), (*ov2)(ii) + (*ov2)(ii));
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_compound_assignment_subtract) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov1(new lorenz95::ObsVec1D(*obstable_));
    ov1->read("ObsVal");
    boost::scoped_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1, true));

    *ov1 -= *ov2;

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov1)(ii), 0.0);
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_compound_assignment_multiply) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov1(new lorenz95::ObsVec1D(*obstable_));
    ov1->read("ObsVal");
    boost::scoped_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1, true));

    *ov1 *= *ov2;

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov1)(ii), (*ov2)(ii) * (*ov2)(ii));
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_compound_assignment_divide) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov1(new lorenz95::ObsVec1D(*obstable_));
    ov1->read("ObsVal");
    boost::scoped_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1, true));

    *ov1 /= *ov2;

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov1)(ii), 1.0);
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_zero) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov(new lorenz95::ObsVec1D(*obstable_));
    ov->read("ObsVal");

    ov->zero();

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov)(ii), 0.0);
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_axpy) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov1(new lorenz95::ObsVec1D(*obstable_));
    ov1->read("ObsVal");
    boost::scoped_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1, true));

    double mult = 2.00;

    ov1->axpy(mult, *ov2);

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov1)(ii), (*ov2)(ii) + mult * (*ov2)(ii));
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_invert) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov1(new lorenz95::ObsVec1D(*obstable_));
    ov1->read("ObsVal");
    boost::scoped_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1, true));

    ov1->invert();

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK_EQUAL((*ov1)(ii), 1.0 / (*ov2)(ii));
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_random) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov(new lorenz95::ObsVec1D(*obstable_));

    ov->random();

    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      BOOST_CHECK((*ov)(ii) != 0);
    }
  }
// -----------------------------------------------------------------------------
  BOOST_AUTO_TEST_CASE(test_ObsVec1D_dot_product_with) {
    boost::scoped_ptr<lorenz95::ObsVec1D> ov1(new lorenz95::ObsVec1D(*obstable_));
    ov1->read("ObsVal");
    boost::scoped_ptr<lorenz95::ObsVec1D> ov2(new lorenz95::ObsVec1D(*ov1, true));

    double result = ov1->dot_product_with(*ov2);

    double check = 0.0;
    for (unsigned int ii = 0; ii < obstable_->nobs(); ++ii) {
      check += (*ov2)(ii) * (*ov2)(ii);
    }

    BOOST_CHECK_CLOSE(result, check, 1.0e-8);
  }
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END()
// -----------------------------------------------------------------------------

}  // namespace test
