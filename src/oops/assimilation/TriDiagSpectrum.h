/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_TRIDIAGSPECTRUM_H_
#define OOPS_ASSIMILATION_TRIDIAGSPECTRUM_H_

#include <vector>

extern "C" {
  void FtnTriDiagSpectrum(const int &, const double *, const double *,
                          double *, double *);
}

namespace oops {

void TriDiagSpectrum(const std::vector<double> & diag, const std::vector<double> & sub,
                     std::vector<double> & evals,
                     std::vector< std::vector<double> > & evecs) {
  const unsigned nn = diag.size();
  ASSERT(sub.size() == nn);
  evals.resize(nn);
  evecs.resize(nn);

  if (nn > 0) {
    std::vector<double> ftnev(nn*nn);

    FtnTriDiagSpectrum(nn, &diag[0], &sub[0], &evals[0], &ftnev[0]);

    unsigned ii = 0;
    for (unsigned jv = 0; jv < nn; ++jv) {
      evecs[jv].resize(nn);
      for (unsigned jj = 0; jj < nn; ++jj) {
        evecs[jv][jj] = ftnev[ii];
        ++ii;
      }
    }
  }
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_TRIDIAGSPECTRUM_H_
