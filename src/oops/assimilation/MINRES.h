/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_MINRES_H_
#define OOPS_ASSIMILATION_MINRES_H_

#include <cmath>
#include <iostream>
#include <vector>

#include "util/dot_product.h"
#include "util/formats.h"
#include "util/Logger.h"

namespace oops {

/*! \file MINRES.h
 * \brief MINRES solver for Ax=b.
 *
 * MINRES solver for Ax=b.(based on implementation following
 * C. C. Paige and M. A. Saunders, 1975)

 * A must be square and symmetric. A can be indefinite.
 * A symmetric positive-definite preconditioner
 * must be supplied that, given a
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
double MINRES(VECTOR & x, const VECTOR & b,
              const AMATRIX & A, const PMATRIX & precond,
              const int maxiter, const double tolerance) {
  VECTOR r(x);
  VECTOR work(x);
  VECTOR v(x);

  r = b;
  double xnrm2 = dot_product(x, x);
  if (xnrm2 != 0) {
    A.multiply(x, work);  // sx = Ax
    r -= work;  // r = b - Ax
  }
  work *= 0;
  VECTOR work2(work);
  VECTOR work1(work);
  VECTOR r1(r);
  VECTOR r2(r);
  VECTOR y(r);

  precond.multiply(r, y);

  double normReduction = 1.0;
  double ynrm2 = sqrt(dot_product(y, y));
  double beta  = sqrt(dot_product(y, r));

  double oldb = 0;
  double epsln = 0;
  double cs = -1.0;
  double sn = 0;
  double dbar = 0;
  double phibar = beta;

  int jiter;
  // MINRES iteration
  Log::info() << std::endl;
  for (jiter = 0; jiter < maxiter; ++jiter) {
    Log::info() << " MINRES Starting Iteration " << jiter+1 << std::endl;

    v = y;
    v *= 1/beta;
    A.multiply(v, y);

    if (jiter > 0) {
      y.axpy(-beta/oldb, r1);
    }

    double alpha = dot_product(y, v);
    y.axpy(-alpha/beta, r2);
    r1 = r2;
    r2 = y;

    precond.multiply(r2, y);

    oldb = beta;
    beta  = sqrt(dot_product(y, r2));

    // Apply Givens Rotation
    double oldeps = epsln;
    double delta = cs*dbar + sn*alpha;
    double gbar = sn*dbar - cs*alpha;
    epsln = sn*beta;
    dbar = -cs*beta;

    // Compute Givens rotation matrix parameters
    double gamma = sqrt(gbar*gbar + beta*beta);  // update this parameter

    cs = gbar/gamma;
    sn = beta/gamma;
    double phi = cs*phibar;
    phibar = sn*phibar;

    // Update x
    double denom = 1/gamma;
    work1 = work2;
    work2 = work;
    work = v;
    work *= denom;
    work.axpy(-oldeps*denom, work1);
    work.axpy(-delta*denom, work2);
    x.axpy(phi, work);

    normReduction = phibar/ynrm2;

    Log::info() << "MINRES end of iteration " << jiter+1 << ". PNorm reduction= "
                << util::full_precision(normReduction) << std::endl << std::endl;

    if (normReduction <= tolerance) {
        Log::info() << "MINRES: Achieved required reduction in residual norm." << std::endl;
        jiter += 1;
        break;
    }
  }

  Log::info() << "MINRES: end" << std::endl;

  return normReduction;
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_MINRES_H_
