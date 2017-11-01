/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_BASE_TESTSUITEOPOBS_H_
#define TEST_BASE_TESTSUITEOPOBS_H_

#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "util/Duration.h"
#include "util/dot_product.h"

namespace test {

// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(test_h_obs_type) {
  BOOST_TEST_MESSAGE("[MSG] Test observation operator for each type");
  const util::Config cfConf(getConfig(), "cost_function");
  util::Duration fclen(cfConf.getData("window_length"));
  util::DateTime winbgn(cfConf.getData("window_begin"));
  util::DateTime winend(winbgn);
  winend += fclen;

  const ControlVariable_ xx(J().jb().getBackground());
  const util::Config biasConfig(getConfig(), "ObsBias", true);
  ObsAuxCtrl_ ybias(biasConfig);

  unsigned nobs = cfConf.getElementSize("Jo");
  for (unsigned jobs = 0; jobs < nobs; ++jobs) {
    ModelState_ xstart(xx.state()[0], getModelConfig());
    const util::Config joConfig(cfConf, "Jo", jobs);
    const util::Config obsConfig(joConfig, "Observation");
    Observations_ yy(obsConfig, winbgn, winend);
    boost::scoped_ptr<Departures_> ref(yy.newDepartures("ombg"));

    boost::shared_ptr<Observer_> pp(new Observer_(yy, ybias));
    oops::PostProcessor<State_> post;
    post.enrollProcessor(pp);

    xstart.forecast(getModelErr(), fclen, post);
    boost::shared_ptr<Observations_> yequ(pp->release());

    Departures_ dep = yy - (*yequ);
    dep.save("oman");

    dep -= *ref;
    double dist = dot_product(dep, dep);

    std::ostringstream oss;
    oss << std::setprecision(20) << "[MSG]   Obs type  (" << jobs << " nbObs= " << dep.numberOfObs()
        << " ) dist[ H(OOPS), H(SCREENING) ] = " << dist;
    BOOST_TEST_MESSAGE(oss.str());
    BOOST_CHECK_SMALL(dist, tolAD());
  }
}

}  // namespace test

#endif  // TEST_BASE_TESTSUITEOPOBS_H_
