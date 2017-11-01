/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_PCG_H_
#define OOPS_ASSIMILATION_PCG_H_

#include <cmath>
#include <vector>

#include "util/dot_product.h"
#include "util/formats.h"
#include "util/Logger.h"

namespace oops {

/*! \file PCG.h
 * \brief Preconditioned Conjugate Gradients solver.
 *
 * This solver is based on the standard Preconditioned Conjugate
 * Gradients solver for Ax=b
 *
 * A must be square, symmetric, positive definite.
 * A preconditioner must be supplied that, given a vector q, returns an
 * approximate solution of Ap=q.
 *
 * On entry:
 * -    x       =  starting point, \f$ x_0 \f$.
 * -    b       = right hand side.
 * -    A       = \f$ A \f$.
 * -    precond = preconditioner \f$ P \approx (A)^{-1} \f$.
 *
 * On exit, x will contain the solution \f$ x \f$

 *  The return value is the achieved reduction in preconditioned residual norm.
 *
 *  Iteration will stop if the maximum iteration limit "maxiter" is reached
 *  or if the preconditioned residual norm reduces by a factor of "tolerance".
 *
 *  VECTOR must implement:
 *  - dot_product
 *  - operator(=)
 *  - operator(+=),
 *  - operator(-=)
 *  - operator(*=) [double * VECTOR],
 *  - axpy
 *
 *  Each of AMATRIX and PMATRIX must implement a method:
 *  - void multiply(const VECTOR&, VECTOR&) const
 *
 *  which applies the matrix to the first argument, and returns the
 *  matrix-vector product in the second. (Note: the const is optional, but
 *  recommended.)
 */

template <typename VECTOR, typename AMATRIX, typename PMATRIX>
double PCG(VECTOR & x, const VECTOR & b,
            const AMATRIX & A, const PMATRIX & precond,
            const int maxiter, const double tolerance ) {
  VECTOR ap(x);
  VECTOR p(x);
  VECTOR r(x);
  VECTOR s(x);
  VECTOR v(x);
  VECTOR z(x);

  std::vector<VECTOR> vVEC;  // required for re-orthogonalization
  std::vector<VECTOR> zVEC;  // required for re-orthogonalization

  // Initial residual r = b - Ax
  r = b;
  double xnrm2 = dot_product(x, x);
  if (xnrm2 != 0) {
    A.multiply(x, s);
    r -= s;
  }

  // s = precond r
  precond.multiply(r, s);

  double dotRr0  = dot_product(r, s);
  double normReduction = 1.0;
  double rdots = dotRr0;
  double rdots_old = 0.0;

  v = r;
  v  *= 1/sqrt(dotRr0);
  z = s;
  z  *= 1/sqrt(dotRr0);

  vVEC.push_back(v);
  zVEC.push_back(z);

  Log::info() << std::endl;
  for (int jiter = 0; jiter < maxiter; ++jiter) {
    Log::info() << " PCG Starting Iteration " << jiter+1 << std::endl;

    if (jiter == 0) {
      p  = s;
    } else {
      double beta = dot_product(s, r)/rdots_old;
      Log::debug() << "PCG beta = " << beta << std::endl;

      p *= beta;
      p += s;      // p = s + beta*p
    }

    A.multiply(p, ap);  // ap = A*p

    double alpha = rdots/dot_product(p, ap);

    x.axpy(alpha, p);    // x = x + alpha*p
    r.axpy(-alpha, ap);  // r = r - alpha*ap

    // Re-orthogonalization
    for (int iiter = 0; iiter < jiter; ++iiter) {
      double proj = dot_product(r, zVEC[iiter]);
      r.axpy(-proj, vVEC[iiter]);
    }

    precond.multiply(r, s);

    rdots_old = rdots;
    rdots = dot_product(r, s);

    v = r;
    v  *= 1/sqrt(rdots);
    z = s;
    z  *= 1/sqrt(rdots);
    vVEC.push_back(v);
    zVEC.push_back(z);

    normReduction = sqrt(rdots/dotRr0);

    Log::info() << "PCG end of iteration " << jiter+1 << ". Norm reduction= "
                << util::full_precision(normReduction) << std::endl << std::endl;

    if (normReduction < tolerance) {
      Log::info() << "PCG: Achieved required reduction in residual norm." << std::endl;
      break;
    }
  }

  Log::info() << "PCG: end" << std::endl;

  return normReduction;
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_PCG_H_
