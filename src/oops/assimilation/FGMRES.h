/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_FGMRES_H_
#define OOPS_ASSIMILATION_FGMRES_H_

#include <cmath>
#include <iostream>
#include <vector>

#include "oops/assimilation/rotmat.h"
#include "oops/assimilation/UpTriSolve.h"
#include "util/dot_product.h"
#include "util/formats.h"
#include "util/Logger.h"

namespace oops {

/*! \file FGMRES.h
 * \brief FGMRES solver for Ax=b.
 *
 * Flexible GMRES solver for Ax=b.(based on implementation following
 * Youcef Saad,SIAM Journal on Scientific Computing, Volume 14 Issue 2,
 * March 1993)

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
double FGMRES(VECTOR & x, const VECTOR & b,
              const AMATRIX & A, const PMATRIX & precond,
              const int maxiter, const double tolerance) {
  std::vector<VECTOR> V;
  std::vector<VECTOR> Z;
  std::vector< std::vector<double> > H;
  std::vector<double> cs(maxiter+1, 0);
  std::vector<double> sn(maxiter+1, 0);
  std::vector<double> s;
  std::vector<double> y(maxiter+1, 0);

  VECTOR r(x);
  VECTOR z(x);
  VECTOR work(x);
  VECTOR work2(x);
  VECTOR w(x);
  VECTOR xinit(x);

  r = b;
  double xnrm2 = dot_product(x, x);
  if (xnrm2 != 0) {
    A.multiply(x, work);
    r -= work;  // r = b - Ax
  }

  double bnrm2  = sqrt(dot_product(b, b));
  double rnrm2  = sqrt(dot_product(r, r));
  double normReduction = rnrm2 / bnrm2;

  // Initialiaze (maxiter + 1) by maxiter matrix H
  H.resize(maxiter);
  for (int ii = 0; ii <= maxiter-1; ii++) {
    H[ii].resize(maxiter + 1);
    for (int jj = 0; jj <= maxiter; jj++) {
       H[ii][jj] = 0;
    }
  }
  work = r;
  work *= 1/rnrm2;
  V.push_back(work);
  s.push_back(rnrm2);

  int jiter;
  // FGMRES iteration
  Log::info() << std::endl;
  for (jiter = 0; jiter < maxiter; ++jiter) {
    Log::info() << " FGMRES Starting Iteration " << jiter+1 << std::endl;

    precond.multiply(V[jiter], z);
    Z.push_back(z);
    A.multiply(z, w);

    double avnrm2 = sqrt(dot_product(w, w));

    // Modified Gram-Schmidt
    for (int jj = 0; jj <= jiter; ++jj) {
      H[jiter][jj] = dot_product(V[jj], w);
      w.axpy(-H[jiter][jj], V[jj]);
    }
    H[jiter][jiter+1] = sqrt(dot_product(w, w));
    double av2nrm2 = H[jiter][jiter+1];

    // Re-orthogonalize if necessary
    if (avnrm2 + 0.001*av2nrm2 == avnrm2) {
      for (int jj = 0; jj <= jiter; ++jj) {
        double hr = dot_product(V[jj], w);
        H[jiter][jj] += hr;
        w.axpy(-hr, V[jj]);
      }
      H[jiter][jiter+1] = sqrt(dot_product(w, w));
    }

    // Check breakdown
    if (H[jiter][jiter+1] != 0) {
      w *= 1/H[jiter][jiter+1];
    }
    V.push_back(w);

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

    H[jiter][jiter] = cs[jiter]*H[jiter][jiter] + sn[jiter]*H[jiter][jiter+1];
    H[jiter][jiter+1] = 0.0;

    double temp = cs[jiter]*s[jiter];
    s.push_back(-sn[jiter]*s[jiter]);
    s[jiter] = temp;

    // Calculate the solution
    UpTriSolve(H, s, y, jiter);  // H s = y
    x = xinit;
    for (int jj = 0; jj < jiter+1; ++jj) {
      x.axpy(y[jj], Z[jj]);
    }

    normReduction = std::abs(s[jiter+1])/bnrm2;
    Log::info() << "FGMRES end of iteration " << jiter+1 << ". PNorm reduction= "
                << util::full_precision(normReduction) << std::endl << std::endl;

    if (normReduction <= tolerance) {
        Log::info() << "FGMRES: Achieved required reduction in residual norm." << std::endl;
        jiter += 1;
        break;
    }
  }

  // Calculate the solution
  UpTriSolve(H, s, y, jiter);  // H s = y
  x = xinit;
  for (int jj = 0; jj < jiter; ++jj) {
    x.axpy(y[jj], Z[jj]);
  }

  Log::info() << "FGMRES: end" << std::endl;

  return normReduction;
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_FGMRES_H_
