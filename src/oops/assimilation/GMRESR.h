/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_GMRESR_H_
#define OOPS_ASSIMILATION_GMRESR_H_

#include <cmath>
#include <vector>

#include "util/dot_product.h"
#include "util/formats.h"
#include "util/Logger.h"

namespace oops {

/*! \file GMRESR.h
 * \brief GMRESR solver for Ax=b.
 *
 * GMRESR solver for Ax=b. (H.A. Van der Vorst and C. Vuik, 1994,
 * Numerical Linear Algebra with Applications, 1(4), 369-386.)
 * A must be square, but need not be symmetric. (For a symmetric matrix,
 * and constant preconditioner, GMRESR is simply PCG with full
 * orthogonalisation.) A preconditioner must be supplied that, given a
 * vector q, returns an approximate solution of Ap=q. The preconditioner
 * can be variable.
 *
 * On entry:
 * -    x       =  starting point, \f$ X_0 \f$.
 * -    b       = right hand side.
 * -    A       = \f$ A \f$.
 * -    precond = preconditioner \f$ F_k \approx (A)^{-1} \f$.
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
double GMRESR(VECTOR & xx, const VECTOR & bb,
              const AMATRIX & A, const PMATRIX & precond,
              const int maxiter, const double tolerance) {
  const double stagthresh = 1.0e-3;
  const double smallres = 1.0e-6;

  std::vector<VECTOR> cvecs;
  std::vector<VECTOR> zvecs;
  VECTOR cc(xx);
  VECTOR zz(xx);

  VECTOR rr(bb);
  A.multiply(xx, zz);  // zz=Axx
  rr -= zz;

  double normReduction = 1.0;
  double rrnorm = sqrt(dot_product(rr, rr));
  double cdotr = rrnorm;
  const double rrnorm0 = rrnorm;

  if (rrnorm > smallres) {
    Log::info() << std::endl;
    for (int jiter = 0; jiter < maxiter; ++jiter) {
      Log::info() << " GMRESR Starting Iteration " << jiter+1 << std::endl;

      // test for stagnation
      if (std::abs(cdotr) < stagthresh*rrnorm) {
        Log::info() << "GMRESR stagnated. Doing an LSQR step." << std::endl;
        A.multiply(rr, zz);  // should be A^T, but we assume A is symmetric.
      } else {
        precond.multiply(rr, zz);  // returns zz as approximate solution of Azz=rr
      }

      A.multiply(zz, cc);  // cc=Azz

      for (int jj = 0; jj < jiter; ++jj) {
        double alpha = -dot_product(cvecs[jj], cc);
        cc.axpy(alpha, cvecs[jj]);  // cc = cc - <c[jj],cc>*c[jj];
        zz.axpy(alpha, zvecs[jj]);  // zz = zz - <c[jj],cc>*z[jj];
      }

      double ccnorm = sqrt(dot_product(cc, cc));
      cvecs.push_back(cc);         // c[jiter]=cc
      cvecs[jiter] *= 1.0/ccnorm;
      zvecs.push_back(zz);         // z[jiter]=zz
      zvecs[jiter] *= 1.0/ccnorm;

      cdotr = dot_product(cvecs[jiter], rr);
      rr.axpy(-cdotr, cvecs[jiter]);
      xx.axpy(cdotr, zvecs[jiter]);

      rrnorm = sqrt(dot_product(rr, rr));
      normReduction = rrnorm/rrnorm0;
      Log::info() << "GMRESR end of iteration " << jiter+1 << ". Norm reduction= "
                  << util::full_precision(normReduction) << std::endl << std::endl;

      if (normReduction < tolerance) {
        Log::info() << "GMRESR: Achieved required reduction in residual norm." << std::endl;
        break;
      }
    }
  }

  return normReduction;
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_GMRESR_H_
