/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_DRPCGMINIMIZER_H_
#define OOPS_ASSIMILATION_DRPCGMINIMIZER_H_

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/BMatrix.h"
#include "oops/assimilation/CMatrix.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DRMinimizer.h"
#include "oops/assimilation/HtRinvHMatrix.h"
#include "oops/assimilation/MinimizerUtils.h"
#include "oops/assimilation/QNewtonLMP.h"
#include "oops/util/dot_product.h"
#include "oops/util/formats.h"
#include "oops/util/Logger.h"

namespace oops {

/// DRPCG Minimizer
/*!
 * \brief Derber-Rosati Preconditioned Conjugate Gradients solver.
 *
 * This solver is based on the standard Preconditioned Conjugate
 * Gradients solver for Ax=b (G. H. Golub and C. F. Van Loan, Matrix
 * Computations), and on the Derber and Rosati double
 * PCG algorithm (J. Derber and A. Rosati, 1989, J. Phys. Oceanog. 1333-1347).
 * For details see S. Gurol, PhD Manuscript, 2013.
 * It solves \f$ Ax=b\f$ for the particular case \f$ A=B^{-1}+C\f$,
 * without requiring the application of \f$ B^{-1}\f$. This algorithm
 * is similar to DRIPCG except it includes standard PCG instead IPCG
 * and stopping criteria is based on the preconditioner norm.
 *
 * A must be square, symmetric, positive definite.
 *
 * A preconditioner must be supplied that, given a vector q, returns an
 * approximation to \f$ (AB)^{-1} q\f$. Possible preconditioning
 * is detailed in S. Gurol, PhD Manuscript, 2013.
 * Note that the traditional \f$ B\f$-preconditioning corresponds to
 * precond=\f$I\f$.
 *
 * On entry:
 * -    dx      =  starting point, \f$ dx_{0} \f$.
 * -    dxh     = \f$ B^{-1} dx_{0} \f$.
 * -    rr      = \f$ (sum B^-1 dx^{b}_{i} + ) H^T R^{-1} d \f$
 * -    B       = \f$ B \f$.
 * -    C       = \f$ C \f$.
 * -    precond = preconditioner \f$ F_k \approx (AB)^{-1} \f$.
 *
 * On exit, dx will contain the solution \f$ dx \f$ and dxh will contain
 * \f$ B^{-1} dx\f$.
 *  The return value is the achieved reduction in preconditioned residual norm.
 *
 *  Iteration will stop if the maximum iteration limit "maxiter" is reached
 *  or if the residual norm reduces by a factor of "tolerance".
 *
 *  Each matrix must implement a method:
 *  - void multiply(const VECTOR&, VECTOR&) const
 *
 *  which applies the matrix to the first argument, and returns the
 *  matrix-vector product in the second. (Note: the const is optional, but
 *  recommended.)
 */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class DRPCGMinimizer : public DRMinimizer<MODEL, OBS> {
  typedef BMatrix<MODEL, OBS>             Bmat_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef HtRinvHMatrix<MODEL, OBS>       HtRinvH_;
  typedef CMatrix<MODEL, OBS>             Cmat_;

 public:
  const std::string classname() const override {return "DRPCGMinimizer";}
  DRPCGMinimizer(const eckit::Configuration &, const CostFct_ &);
  ~DRPCGMinimizer() {}

 private:
  double solve(CtrlInc_ &, CtrlInc_ &, CtrlInc_ &, const Bmat_ &, const HtRinvH_ &,
               const double, const double, const int, const double) override;
  QNewtonLMP<CtrlInc_, Bmat_, Cmat_> lmp_;
};

// =============================================================================

template<typename MODEL, typename OBS>
DRPCGMinimizer<MODEL, OBS>::DRPCGMinimizer(const eckit::Configuration & conf, const CostFct_ & J)
  : DRMinimizer<MODEL, OBS>(J), lmp_(conf)
{}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double DRPCGMinimizer<MODEL, OBS>::solve(CtrlInc_ & dx, CtrlInc_ & dxh, CtrlInc_ & rr,
                                   const Bmat_ & B, const HtRinvH_ & HtRinvH,
                                   const double costJ0Jb, const double costJ0JoJc,
                                   const int maxiter, const double tolerance) {
  // dx   increment
  // dxh  B^{-1} dx
  // rr   (sum B^{-1} dx_i^{b} +) G^T H^{-1} d

  CtrlInc_ qq(dxh);
  CtrlInc_ pp(dxh);
  CtrlInc_ hh(dxh);
  CtrlInc_ zz(dxh);
  CtrlInc_ pr(dxh);
  CtrlInc_ r0(dxh);
  CtrlInc_ ww(dxh);

  // vectors for re-orthogonalization
  std::vector<CtrlInc_> vvecs;
  std::vector<CtrlInc_> zvecs;
  std::vector<double> scals;
  // reserve space in vectors to avoid extra copies
  vvecs.reserve(maxiter+1);
  zvecs.reserve(maxiter+1);
  scals.reserve(maxiter+1);

  // J0
  const double costJ0 = costJ0Jb + costJ0JoJc;

  // r_{0}
  r0 = rr;

  // r_{0}^T r_{0}
  Log::info() << "normr0 " << dot_product(rr, rr) << std::endl;

  // Set ObsBias part of the preconditioner
  lmp_.updateObsBias(std::make_unique<Cmat_>(B.obsAuxCovariance()));

  // z_{0} = B LMP r_{0}
  lmp_.multiply(rr, pr);
  B.multiply(pr, zz);
  // p_{0} = z_{0}
  pp = zz;
  // h_{0} = LMP r_{0}
  hh = pr;

  // r_{0}^T z_{0}
  double dotRr0  = dot_product(rr, zz);
  double normReduction = 1.0;
  double rdots = dotRr0;
  double rdots_old = 0.0;

  vvecs.push_back(rr);
  zvecs.push_back(zz);
  scals.push_back(1.0/dotRr0);

  Log::info() << std::endl;
  for (int jiter = 0; jiter < maxiter; ++jiter) {
    Log::info() << " DRPCG Starting Iteration " << jiter+1 << std::endl;

    if (jiter > 0) {
      // beta_{i} = r_{i+1}^T z_{i+1} / r_{i}^T z_{i}
      double beta = rdots/rdots_old;

      // p_{i+1} = z_{i+1} + beta*p_{i}
      pp *= beta;
      pp += zz;

      // h_{i+1} = LMP r_{i+1} + beta*h_{i}
      hh *= beta;
      hh += pr;
    }

    // q_{i} = h_{i} + H^T R^{-1} H p_{i}
    HtRinvH.multiply(pp, qq);
    qq += hh;

    // alpha_{i} = r_{i}^T z_{i} / q_{i}^T p_{i}
    double rho = dot_product(pp, qq);
    double alpha = rdots/rho;

    // dx_{i+1} = dx_{i} + alpha * p_{i}
    dx.axpy(alpha, pp);
    // dxh_{i+1} = dxh_{i} + alpha * h_{i} ! for diagnosing Jb
    dxh.axpy(alpha, hh);
    // r_{i+1} = r_{i} - alpha * q_{i}
    rr.axpy(-alpha, qq);

    // Compute the quadratic cost function
    // J[dx_{i}] = J[0] - 0.5 dx_{i}^T r_{0}
    double costJ = costJ0 - 0.5 * dot_product(dx, r0);
    // Jb[dx_{i}] = 0.5 dx_{i}^T f_{i}
    double costJb = costJ0Jb + 0.5 * dot_product(dx, dxh);
    // Jo[dx_{i}] + Jc[dx_{i}] = J[dx_{i}] - Jb[dx_{i}]
    double costJoJc = costJ - costJb;

    // Re-orthogonalization
    for (int jj = 0; jj < jiter; ++jj) {
      double proj = scals[jj] * dot_product(rr, zvecs[jj]);
      rr.axpy(-proj, vvecs[jj]);
    }

    // z_{i+1} = B LMP r_{i+1}
    lmp_.multiply(rr, pr);
    B.multiply(pr, zz);

    // r_{i}^T z_{i}
    rdots_old = rdots;
    // r_{i+1}^T z_{i+1}
    rdots = dot_product(rr, zz);

    Log::info() << "rdots " << rdots << "      iteration    " << jiter << std::endl;

    // r_{i+1}^T z_{i+1} / r_{0}^T z_{0}
    normReduction = sqrt(rdots/dotRr0);

    Log::info() << "DRPCG end of iteration " << jiter+1 << std::endl;
    printNormReduction(jiter+1, sqrt(rdots), normReduction);
    printQuadraticCostFunction(jiter+1, costJ, costJb, costJoJc);

    // Save the pairs for preconditioning
    lmp_.push(pp, hh, qq, rho);

    if (normReduction < tolerance) {
      Log::info() << "DRPCG: Achieved required reduction in residual norm." << std::endl;
      break;
    }

    vvecs.push_back(rr);
    zvecs.push_back(zz);
    scals.push_back(1.0/rdots);
  }

// Generate the (second-level) Limited Memory Preconditioner
  lmp_.update(B);

  return normReduction;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_DRPCGMINIMIZER_H_
