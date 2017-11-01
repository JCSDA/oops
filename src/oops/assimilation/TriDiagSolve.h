/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_TRIDIAGSOLVE_H_
#define OOPS_ASSIMILATION_TRIDIAGSOLVE_H_

#include <vector>

namespace oops {

void TriDiagSolve(const std::vector<double> & diag, const std::vector<double> & sub,
                  const std::vector<double> & rhs, std::vector<double> & sol) {
  const double n = rhs.size();
  ASSERT(sub.size() == n-1);
  ASSERT(diag.size() == n);
  sol.resize(n);
  std::vector<double> c(sub);
  std::vector<double> d(rhs);

  // Forward substitution
  c[0] = c[0] / diag[0];
  d[0] = d[0] / diag[0];
  if (n > 2) {
    for (int iiter = 1; iiter <= n-2; ++iiter) {
      double denom = diag[iiter]-sub[iiter-1]*c[iiter-1];
      c[iiter] = c[iiter]/denom;
      d[iiter] = (d[iiter]-sub[iiter-1]*d[iiter-1])/denom;
    }
  }
  d[n-1] = (d[n-1]-sub[n-2]*d[n-2])/(diag[n-1]-sub[n-2]*c[n-2]);
  // Backward substitution
  sol[n-1] = d[n-1];
  for (int iiter = n-2; iiter >= 0; --iiter) {
    sol[iiter] = (d[iiter]-(c[iiter]*sol[iiter+1]));
  }
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_TRIDIAGSOLVE_H_
