/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_BASE_TESTSUITEOPOBSTL_H_
#define TEST_BASE_TESTSUITEOPOBSTL_H_

#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace test {

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_tl_obs_type) {
  BOOST_TEST_MESSAGE("[MSG] Test Tangent linear observation operator for each type");

  const ControlVariable_ xx(J().jb().getBackground());
  boost::scoped_ptr<CtrlInc_> dx(new CtrlInc_(J().jb()));
  randomize(*dx);
  (*dx) *= deltaX();

  for (unsigned jj = 0; jj < J().nterms(); ++jj) {
    boost::shared_ptr<Observations_> H_x = multiply_Hnl(xx, jj);
    boost::shared_ptr<Departures_> H_dx = multiply_Htl(*dx, jj);

    double eps = 10.0;
    for (int lambda = 1; lambda <= 10; lambda++) {
      eps /= 10.0;

      boost::scoped_ptr<Departures_> eps_H_dx(new Departures_(*H_dx));
      (*eps_H_dx) *= eps;

      boost::scoped_ptr<CtrlInc_> eps_dx(new CtrlInc_(*dx));
      *eps_dx *= eps;

      ControlVariable_ x_plus_eps_dx(xx);
      J().addIncrement(x_plus_eps_dx, *eps_dx);
      boost::shared_ptr<Observations_> H_x_plus_eps_dx = multiply_Hnl(x_plus_eps_dx, jj);

      Departures_ ratio = (*H_x_plus_eps_dx) - (*H_x);
      ratio /= (*eps_H_dx);

      /* This test makes sense only with a single obs file */
      double scal_ratio = sqrt(dot_product(ratio, ratio));

      std::ostringstream oss;
      oss << std::setprecision(20) << "[MSG]   Obs type  (" << jj << ") epsilon=" << eps
          << " ( Hnl(x+eps.dx) - Hnl(x) ) / eps.Htl(dx) = " << scal_ratio;
      BOOST_TEST_MESSAGE(oss.str());

//    BOOST_CHECK_SMALL(scal_ratio, tolAD());
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_BASE_TESTSUITEOPOBSTL_H_
