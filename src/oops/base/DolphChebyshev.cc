/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/base/DolphChebyshev.h"

#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

DolphChebyshev::DolphChebyshev(const eckit::Configuration & config) {
  tau_ = util::Duration(config.getString("cutoff"));
}

// -----------------------------------------------------------------------------

std::map< util::DateTime, double > DolphChebyshev::setWeights(const util::DateTime & bgn,
                                                              const util::DateTime & end,
                                                              const util::Duration & dt) {
  const double pi = 4.0*std::atan(1.0);
  const util::Duration window(end-bgn);
  const int nstep = window.toSeconds() / dt.toSeconds();
  ASSERT(window.toSeconds() == dt.toSeconds()*nstep);

  const int M = nstep/2;
  const int N = 2*M + 1;

  const int tt = tau_.toSeconds() / dt.toSeconds();
  ASSERT(tau_.toSeconds() == dt.toSeconds()*tt);
  ASSERT(tt > 1);
  const double thetas = 2.0 * pi / tt;
  const double x0 = 1.0 / std::cos(thetas/2.0);
  const double rr = 1.0 / std::cosh(nstep * std::acosh(x0));

  std::vector<double> w(M+1);
  for (int n = 0; n <= M; ++n) {
    double tn = 2.0 * pi * n / N;
    double sum = 0.0;
    for (int m = 1; m <= M; ++m) {
      double xx = x0 * std::cos(pi * m / N);
      double t2m = 0.0;
      double tnm2 = 1.0;
      double tnm1 = xx;
      for (int kk = 2; kk <= 2 * M; ++kk) {
        t2m = 2.0 * xx * tnm1 - tnm2;
        tnm2 = tnm1;
        tnm1 = t2m;
      }
      sum += t2m * std::cos(m*tn);
    }
    w[n] = (1.0 + 2.0 * rr * sum) / N;
  }

  std::map< util::DateTime, double > weights;
  double checksum = 0.0;
  util::DateTime now(bgn);
  for (int jj = -M; jj <= M; ++jj) {
    int n = std::abs(jj);
    weights[now] = w[n];
    checksum += w[n];
    now += dt;
  }
  ASSERT(now == end+dt);
  ASSERT(std::abs(checksum-1.0) < 1.0e-8);

  return weights;
}

// -----------------------------------------------------------------------------

}  // namespace oops

