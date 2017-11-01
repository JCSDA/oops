/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_UPTRISOLVE_H_
#define OOPS_ASSIMILATION_UPTRISOLVE_H_

#include <vector>

namespace oops {

/*! \file UpTriSolve.h
 * \brief Solves the upper triangular system sol = H \ rhs
*/

void UpTriSolve(const std::vector <std::vector<double> > & H,
                const std::vector <double> & rhs, std::vector <double> & sol,
                const int & dim) {
  // Backward solve
  for (int iiter = dim - 1; iiter >= 0; iiter--) {
    sol[iiter] = rhs[iiter];
    for (int jiter = iiter + 1; jiter < dim; jiter++) {
      sol[iiter] -= H[jiter][iiter] * sol[jiter];
    }
    sol[iiter] *= 1/H[iiter][iiter];
  }
}
}  // namespace oops

#endif  // OOPS_ASSIMILATION_UPTRISOLVE_H_
