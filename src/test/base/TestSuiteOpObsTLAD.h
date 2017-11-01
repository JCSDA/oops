/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_BASE_TESTSUITEOPOBSTLAD_H_
#define TEST_BASE_TESTSUITEOPOBSTLAD_H_

#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "util/Logger.h"
#include "util/Duration.h"
#include "util/dot_product.h"

namespace test {

// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_adjoint_HMatrix_HtMatrix) {
  BOOST_TEST_MESSAGE(
      "[MSG] Test observation operator adjoint with HMatrix and HtMatrix");

  // Call H_tl
  boost::scoped_ptr<CtrlInc_> dx1(new CtrlInc_(J().jb()));
  randomize(*dx1);

  DepartVector_ dy1;
  dy1.clear();

  const HMatrix_ H(J());
  H.multiply(*dx1, dy1);

  // Call H_ad
  boost::scoped_ptr<CtrlInc_> dx2(new CtrlInc_(J().jb()));
  dx2->zero();

  DepartVector_ dy2;
  randomize(dy2);

  const HtMatrix_ Ht(J());
  Ht.multiply(dy2, *dx2);

  // Check outcomes
  double dx1dx2 = dot_product(*dx1, *dx2);
  double dy1dy2 = dot_product(dy1, dy2);

  std::ostringstream oss;
  oss << std::setprecision(20) << "[MSG]   dx1.dx2 = " << dx1dx2 << " dy1.dy2 = "
  << dy1dy2 << " diff = " << (dx1dx2-dy1dy2);
  BOOST_TEST_MESSAGE(oss.str());

  BOOST_CHECK_CLOSE(dx1dx2, dy1dy2, tolAD());
}

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(test_adjoint_obs_type) {
  BOOST_TEST_MESSAGE(
      "[MSG] Test observation operator adjoint for each type");

  for (unsigned jj = 0; jj < J().nterms(); ++jj) {
    // Call H_tl
    boost::scoped_ptr<CtrlInc_> dx1(new CtrlInc_(J().jb()));
    randomize(*dx1);

    boost::shared_ptr<Departures_> dy1;
    PostProcessorTL_ costtl;
    costtl.enrollProcessor(J().jterm(jj).setupTL(*dx1));

    J().runTLM(*dx1, costtl);
    dy1.reset(dynamic_cast<Departures_*>(costtl.releaseOutputFromTL(0)));

    // Call H_ad
    boost::scoped_ptr<CtrlInc_> dx2(new CtrlInc_(J().jb()));
    J().zeroAD(*dx2);
    boost::shared_ptr<Departures_> dy2(getRandomDepartures(jj));

    PostProcessorAD_ costad;
    costad.enrollProcessor(J().jterm(jj).setupAD(dy2, *dx2));

    J().runADJ(*dx2, costad);

    // Check outcomes
    double dx1dx2 = dot_product(*dx1, *dx2);
    double dy1dy2 = dot_product(*dy1, *dy2);

    std::ostringstream oss;
    oss << std::setprecision(20) << "[MSG]   Obs type  (" << jj << " nbObs= " << dy1->numberOfObs() << " ) : dx1.dx2 = "
    << dx1dx2 << " dy1.dy2 = " << dy1dy2 << " diff = " << (dx1dx2-dy1dy2);
    BOOST_TEST_MESSAGE(oss.str());

    BOOST_CHECK_CLOSE(dx1dx2, dy1dy2, tolAD());
  }
}

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_BASE_TESTSUITEOPOBSTLAD_H_
