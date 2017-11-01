/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_FULLGMRES_H_
#define OOPS_ASSIMILATION_FULLGMRES_H_

#include <cmath>
#include <iostream>
#include <vector>

#include "oops/assimilation/rotmat.h"
#include "oops/assimilation/UpTriSolve.h"
#include "util/dot_product.h"
#include "util/formats.h"
#include "util/Logger.h"

namespace oops {

/*! \file FullGMRES.h
 * \brief FullGMRES solver for Ax=b.
 *
 * GMRES solver for Ax=b.(based on implementation following
 * Saad-Schultz and C. T. Kelley, July 10, 1994)

 * A must be square, but need not be symmetric.
 * A preconditioner must be supplied that, given a
 * vector q, returns an approximate solution of Ap=q.
 *
 * On entry:
 * -    x       =  starting point, \f$ X_0 \f$.
 * -    b       = right hand side.
 * -    A       = \f$ A \f$.
 * -    precond = preconditioner \f$ F_k \approx (A)^{-1} \f$.
 * -    maxiter = maximum number of iterations
 * -    tol     = error tolerance
 *
 * On exit, x will contain the solution.
 * The return value is the achieved reduction in residual norm.
 *
 * Iteration will stop if the maximum iteration limit "maxiter" is reached
 * or if the residual norm reduces by a factor of "tolerance".
 *
 * VECTOR must implement:
 *  - dot_product
 *  - operator(=)
 *  - operator(+=),
 *  - operator(-=)
 *  - operator(*=) [double * VECTOR],
 *  - axpy
 *
 *  AMATRIX and PMATRIX must implement a method:
 *  - void multiply(const VECTOR&, VECTOR&) const
 *
 *  which applies the matrix to the first argument, and returns the
 *  matrix-vector product in the second. (Note: the const is optional, but
 *  recommended.)
 */

template <typename VECTOR, typename AMATRIX, typename PMATRIX>
double FullGMRES(VECTOR & xx, const VECTOR & bb, const AMATRIX & A,
                 const PMATRIX & precond, const int maxiter,
                 const double tolerance, std::vector<VECTOR> & pqVEC,
                 std::vector<VECTOR> & xyVEC) {
  std::vector<VECTOR> VV;
  std::vector< std::vector<double> > H;
  std::vector<double> cs(maxiter+1, 0);
  std::vector<double> sn(maxiter+1, 0);
  std::vector<double> ss;
  std::vector<double> yy(maxiter+1, 0);

  VECTOR ww(xx);
  VECTOR zz(xx);

  VECTOR rr(bb);
  double xnrm2 = dot_product(xx, xx);
  if (xnrm2 != 0) {
    A.multiply(xx, ww);
    rr -= ww;  // r = b - Ax
  }

  precond.multiply(rr, zz);

  double znrm2 = sqrt(dot_product(zz, zz));
  double normReduction = 1.0;

  zz *= 1/znrm2;
  VV.push_back(zz);
  ss.push_back(znrm2);

  // Initialiaze (maxiter + 1) by maxiter matrix H
  H.resize(maxiter);
  for (int ii = 0; ii <= maxiter-1; ii++) {
    H[ii].resize(maxiter + 1);
    for (int jj = 0; jj <= maxiter; jj++) {
       H[ii][jj] = 0;
    }
  }

  pqVEC.clear();
  xyVEC.clear();

  int jiter;
  // FullGMRES iteration
  Log::info() << std::endl;
  for (jiter = 0; jiter < maxiter; ++jiter) {
    Log::info() << " FullGMRES Starting Iteration " << jiter+1 << std::endl;

    A.multiply(VV[jiter], zz);
    precond.multiply(zz, ww);
    if (jiter > 1) {
      xyVEC.push_back(VV[jiter]);
      pqVEC.push_back(zz);
    }

    double avnrm2 = sqrt(dot_product(ww, ww));

    // Modified Gram-Schmidt
    for (int jj = 0; jj <= jiter; ++jj) {
      H[jiter][jj] = dot_product(VV[jj], ww);
      ww.axpy(-H[jiter][jj], VV[jj]);
    }

    H[jiter][jiter+1] = sqrt(dot_product(ww, ww));

    double av2nrm2 = H[jiter][jiter+1];

    // Re-orthogonalize if necessary
    if (avnrm2 + 0.001*av2nrm2 == avnrm2) {
      for (int jj = 0; jj <= jiter; ++jj) {
        double hr = dot_product(VV[jj], ww);
        H[jiter][jj] += hr;
        ww.axpy(-hr, VV[jj]);
      }
      H[jiter][jiter+1] = sqrt(dot_product(ww, ww));
    }

    // Check breakdown
    if (H[jiter][jiter+1] != 0.0) {
      ww *= 1/H[jiter][jiter+1];
    }

    VV.push_back(ww);

    if (jiter > 0) {
      // Apply Givens Rotation
      for (int jj = 0; jj < jiter; ++jj) {
        double temp = cs[jj]*H[jiter][jj] + sn[jj]*H[jiter][jj+1];
        H[jiter][jj+1] = -sn[jj]*H[jiter][jj] + cs[jj]*H[jiter][jj+1];
        H[jiter][jj] = temp;
      }
    }

    // Compute Givens rotation matrix parameters
    rotmat(H[jiter][jiter], H[jiter][jiter+1], cs[jiter], sn[jiter]);

    H[jiter][jiter]   = cs[jiter]*H[jiter][jiter] + sn[jiter]*H[jiter][jiter+1];
    H[jiter][jiter+1] = 0.0;

    double temp = cs[jiter]*ss[jiter];
    ss.push_back(-sn[jiter]*ss[jiter]);
    ss[jiter] = temp;

    normReduction = std::abs(ss[jiter+1])/znrm2;
    Log::info() << "FullGMRES end of iteration " << jiter+1 << ". PNorm reduction= "
                << util::full_precision(normReduction) << std::endl << std::endl;

    if (normReduction <= tolerance) {
      Log::info() << "FullGMRES: Achieved required reduction in presidual norm." << std::endl;
      jiter += 1;
      break;
    }
  }

  // Calculate the solution
  UpTriSolve(H, ss, yy, jiter);
  for (int jj = 0; jj < jiter; ++jj) {
    xx.axpy(yy[jj], VV[jj]);
  }

  return normReduction;
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_FULLGMRES_H_
