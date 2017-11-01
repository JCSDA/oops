/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_IPCG_H_
#define OOPS_ASSIMILATION_IPCG_H_

#include <cmath>
#include <vector>

#include "util/dot_product.h"
#include "util/formats.h"
#include "util/Logger.h"

namespace oops {

/*! \file IPCG.h
 * \brief Inexact-Preconditioned Conjugate Gradients solver.
 *
 * Golub-Ye Inexact-Preconditioned Conjugate Gradients solver for Ax=b.
 * (G.H. Golub and Q. Ye 1999/00, SIAM J. Sci. Comput. 21(4) 1305-1320.)
 *
 * A must be square, symmetric, positive definite.
 * A preconditioner must be supplied that, given a vector q, returns an
 * approximate solution of Ap=q. The preconditioner can be variable.
 *
 * On entry:
 * -    x       =  starting point, \f$ x_0 \f$.
 * -    b       = right hand side.
 * -    A       = \f$ A \f$.
 * -    precond = preconditioner \f$ F_k \approx (A)^{-1} \f$.
 *
 * On exit, x will contain the solution.
 * The return value is the achieved reduction in residual norm.
 *
 *  Iteration will stop if the maximum iteration limit "maxiter" is reached
 *  or if the residual norm reduces by a factor of "tolerance".
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
double IPCG(VECTOR & x, const VECTOR & b,
            const AMATRIX & A, const PMATRIX & precond,
            const int maxiter, const double tolerance ) {
  VECTOR ap(x);
  VECTOR p(x);
  VECTOR r(x);
  VECTOR s(x);
  VECTOR oldr(x);
  VECTOR w(x);
  VECTOR v(x);   // required for re-orthogonalization
  VECTOR z(x);   // required for re-orthogonalization

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

  double dotRr0  = dot_product(r, r);
  double dotSr0  = dot_product(r, s);
  double normReduction = 1.0;
  double rdots_old = dotSr0;
  double rdots = dotSr0;

  v = r;
  v  *= 1/sqrt(dotSr0);
  z = s;
  z  *= 1/sqrt(dotSr0);

  vVEC.push_back(v);
  zVEC.push_back(z);

  Log::info() << std::endl;
  for (int jiter = 0; jiter < maxiter; ++jiter) {
    Log::info() << " IPCG Starting Iteration " << jiter+1 << std::endl;

    if (jiter == 0) {
      p = s;
    } else {
      w = r;
      w -= oldr;  // w=r-oldr
      double beta = dot_product(s, w)/rdots_old;
      p *= beta;
      p += s;    // p = s + beta*p
    }

    A.multiply(p, ap);  // ap = Ap

    oldr = r;

    double alpha = rdots/dot_product(p, ap);

    x.axpy(alpha, p);    // x = x + alpha*p;
    r.axpy(-alpha, ap);  // r = r - alpha*ap;

    // Re-orthogonalization
    for (int iiter = 0; iiter < jiter; ++iiter) {
      double proj = dot_product(r, zVEC[iiter]);
      r.axpy(-proj, vVEC[iiter]);
    }

    precond.multiply(r, s);  // returns s as approximate solve of As=r

    rdots_old = rdots;
    rdots = dot_product(r, s);

    v = r;
    v  *= 1/sqrt(rdots);
    z = s;
    z  *= 1/sqrt(rdots);
    vVEC.push_back(v);
    zVEC.push_back(z);

    normReduction = sqrt(dot_product(r, r)/dotRr0);
    Log::info() << "IPCG end of iteration " << jiter+1 << ". Norm reduction= "
                << util::full_precision(normReduction) << std::endl << std::endl;

    if (normReduction < tolerance) {
      Log::info() << "IPCG: Achieved required reduction in residual norm." << std::endl;
      break;
    }
  }

  Log::info() << "IPCG: end" << std::endl;

  return normReduction;
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_IPCG_H_
