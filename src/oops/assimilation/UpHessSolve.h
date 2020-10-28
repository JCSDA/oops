/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_UPHESSSOLVE_H_
#define OOPS_ASSIMILATION_UPHESSSOLVE_H_

#include <vector>

#include "oops/assimilation/UpTriSolve.h"

namespace oops {

/*! \file UpHessSolev.h
 *  \brief Solves the linear system H sol = rhs where H is an n-by-n
 *  upper Hessenberg matrix.
 *  Finds the LU factorization of H : H = LU where L is an n-by-n
 *  lower unit triangular and U is an n-by-n upper triangular matrix.
 *  First solves L y = rhs then solves U sol = y
*/

void UpHessSolve(std::vector< std::vector<double> > & UpHess, const std::vector<double> & rhs,
                 std::vector<double> & sol) {
  const double n = rhs.size();
  ASSERT(UpHess.size() == n);
  sol.resize(n);
  std::vector<double> v(n);
  std::vector<double> y(n);

  // Step 1: Compute the LU factorization
  // Note that L = I + diag(v(2:n),-1)
  v[0] = 0;
  for (int ii = 0; ii <= n-2; ++ii) {
    v[ii+1] = UpHess[ii][ii+1]/UpHess[ii][ii];
    for (int jj = ii; jj <= n-1; ++jj) {
      UpHess[jj][ii+1] = UpHess[jj][ii+1] - v[ii+1]*UpHess[jj][ii];
    }
  }
  // Extract the upper triangular part (U part)
  for (int ii = 0; ii <= n-1; ++ii) {
    for (int jj = ii+1; jj <= n-1; ++jj) {
        UpHess[ii][jj] = 0;
    }
  }

  // Step 2: Solve the n-by-n unit lower bidiagonal system  L y = rhs
  y[0] = rhs[0];
  for (int ii = 1; ii <= n-1; ++ii) {
    y[ii] = rhs[ii] - v[ii]*y[ii-1];
  }

  // Step 3: Solve the upper triangular system U sol = y
  UpTriSolve(UpHess, y, sol, n);
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_UPHESSSOLVE_H_
