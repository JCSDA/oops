/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_DRPLANCZOSMINIMIZER_H_
#define OOPS_ASSIMILATION_DRPLANCZOSMINIMIZER_H_

#include <Eigen/Dense>

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "oops/assimilation/BMatrix.h"
#include "oops/assimilation/CMatrix.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DRMinimizer.h"
#include "oops/assimilation/HtRinvHMatrix.h"
#include "oops/assimilation/MinimizerUtils.h"
#include "oops/assimilation/SpectralLMP.h"
#include "oops/assimilation/TriDiagSolve.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "oops/util/printRunStats.h"

namespace oops {

/// DRPLanczos Minimizer
/*!
 * \brief Derber-Rosati Preconditioned Lanczos solver.
 *
 * This solver is the Lanczos version of the DRPCG algorithm
 * It solves \f$ Ax=b\f$ for the particular case \f$ A=B^{-1}+C\f$,
 * without requiring the application of \f$ B^{-1}\f$.
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
 * -    dx      = starting point.
 * -    dxh     = starting point, \f$ B^{-1} dx_{0}\f$.
 * -    rr      = residual at starting point.
 * -    B       = \f$ B \f$.
 * -    C       = \f$ C \f$.
 * -    precond = preconditioner \f$ F_k \approx (AB)^{-1} \f$.
 *
 * On exit, dxh will contain \f$ B^{-1} x\f$ where x is the solution.
 * The return value is the achieved reduction in residual norm.
 *
 * Iteration will stop if the maximum iteration limit "maxiter" is reached
 * or if the residual norm reduces by a factor of "tolerance".
 *
 * Each matrix must implement a method:
 * - void multiply(const VECTOR&, VECTOR&) const
 *
 * which applies the matrix to the first argument, and returns the
 * matrix-vector product in the second. (Note: the const is optional, but
 * recommended.)
 */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class DRPLanczosMinimizer : public DRMinimizer<MODEL, OBS> {
  typedef BMatrix<MODEL, OBS>             Bmat_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef HtRinvHMatrix<MODEL, OBS>       HtRinvH_;
  typedef CMatrix<MODEL, OBS>             Cmat_;

 public:
  const std::string classname() const override {return "DRPLanczosMinimizer";}
  DRPLanczosMinimizer(const eckit::Configuration &, const CostFct_ &);
  ~DRPLanczosMinimizer() {}

 private:
  double solve(CtrlInc_ &, CtrlInc_ &, CtrlInc_ &, const Bmat_ &, const HtRinvH_ &,
               const double, const double, const int, const double) override;

  SpectralLMP<CtrlInc_, Cmat_> lmp_;

  std::vector<std::unique_ptr<CtrlInc_>> hvecs_;
  std::vector<std::unique_ptr<CtrlInc_>> vvecs_;
  std::vector<std::unique_ptr<CtrlInc_>> zvecs_;
  std::vector<double> alphas_;
  std::vector<double> betas_;

  // For diagnostics
  eckit::LocalConfiguration diagConf_;
};

// =============================================================================

template<typename MODEL, typename OBS>
DRPLanczosMinimizer<MODEL, OBS>::DRPLanczosMinimizer(const eckit::Configuration & conf,
                                                     const CostFct_ & J)
  : DRMinimizer<MODEL, OBS>(J), lmp_(conf), hvecs_(), vvecs_(), zvecs_(), alphas_(),
    betas_(), diagConf_(conf) {}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double DRPLanczosMinimizer<MODEL, OBS>::solve(CtrlInc_ & dx, CtrlInc_ & dxh, CtrlInc_ & rr,
                                        const Bmat_ & B, const HtRinvH_ & HtRinvH,
                                        const double costJ0Jb, const double costJ0JoJc,
                                        const int maxiter, const double tolerance) {
  util::printRunStats("DRPLanczos start");
  // dx   increment
  // dxh  B^{-1} dx
  // rr   (sum B^{-1} dx_i^{b} +) G^T H^{-1} d

  CtrlInc_ zz(dxh);
  CtrlInc_ pr(dxh);
  CtrlInc_ vv(rr);

  std::vector<double> ss;
  std::vector<double> dd;

  // Set ObsBias part of the preconditioner
  lmp_.updateObsBias(std::make_unique<Cmat_>(B.obsAuxCovariance()));

  // J0
  const double costJ0 = costJ0Jb + costJ0JoJc;

  lmp_.update(vvecs_, hvecs_, zvecs_, alphas_, betas_);

  // z_{0} = B LMP r_{0}
  lmp_.multiply(vv, pr);
  B.multiply(pr, zz);

  // beta_{0} = sqrt( z_{0}^T r_{0} )
  double beta = sqrt(dot_product(zz, vv));
  const double beta0 = beta;

  // v_{1} = r_{0} / beta_{0}
  vv *= 1/beta;
  // pr_{1} = LMP r_{0} / beta_{0}
  pr *= 1/beta;
  // z_{1} = z_{0} / beta_{0}
  zz *= 1/beta;

  // hvecs[0] = pr_{1} --> required for solution
  hvecs_.emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(pr)));
  // zvecs[0] = z_{1} ---> for re-orthogonalization
  zvecs_.emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(zz)));
  // vvecs[0] = v_{1} ---> for re-orthogonalization
  vvecs_.emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(vv)));

  double normReduction = 1.0;

  Log::info() << std::endl;
  for (int jiter = 0; jiter < maxiter; ++jiter) {
    Log::info() << "DRPLanczos Starting Iteration " << jiter+1 << std::endl;
    util::printRunStats("DRPLanczos iteration " + std::to_string(jiter+1));

    // v_{i+1} = ( pr_{i} + H^T R^{-1} H z_{i} ) - beta * v_{i-1}
    HtRinvH.multiply(zz, vv);
    vv += pr;

    if (jiter > 0) {
      vv.axpy(-beta, *vvecs_[jiter-1]);
    }

    // alpha_{i} = v_{i+1}^T z_{i}
    double alpha = dot_product(zz, vv);

    // v_{i+1} = v_{i+1} - alpha_{i} v_{i}
    vv.axpy(-alpha, *vvecs_[jiter]);  // vv = vv - alpha * v_j

    // Re-orthogonalization
    for (int jj = 0; jj < jiter; ++jj) {
      double proj = dot_product(vv, *zvecs_[jj]);
      vv.axpy(-proj, *vvecs_[jj]);
    }

    // z_{i+1} = B LMP v_{i+1}
    lmp_.multiply(vv, pr);
    B.multiply(pr, zz);

    // beta_{i+1} = sqrt( zz_{i+1}^t, vv_{i+1} )
    beta = sqrt(dot_product(zz, vv));

    // v_{i+1} = v_{i+1} / beta_{i+1}
    vv *= 1/beta;
    // pr_{i+1} = pr_{i+1} / beta_{i+1}
    pr *= 1/beta;
    // z_{i+1} = z_{i+1} / beta_{i+1}
    zz *= 1/beta;

    // hvecs[i+1] =pr_{i+1}
    hvecs_.emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(pr)));
    // zvecs[i+1] = z_{i+1}
    zvecs_.emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(zz)));
    // vvecs[i+1] = v_{i+1}
    vvecs_.emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(vv)));

    alphas_.push_back(alpha);

    if (jiter == 0) {
      ss.push_back(beta0/alpha);
      dd.push_back(beta0);
    } else {
      // Solve the tridiagonal system T_{i} s_{i} = beta0 * e_1
      dd.push_back(beta0*dot_product(*zvecs_[0], vv));
      TriDiagSolve(alphas_, betas_, dd, ss);
    }

    betas_.push_back(beta);

    // Compute the quadratic cost function
    // J[du_{i}] = J[0] - 0.5 s_{i}^T Z_{i}^T r_{0}
    // Jb[du_{i}] = 0.5 s_{i}^T V_{i}^T Z_{i} s_{i}
    double costJ = costJ0;

    double costJb = costJ0Jb;
    for (int jj = 0; jj < jiter+1; ++jj) {
      costJ -= 0.5 * ss[jj] * dot_product(*zvecs_[jj], rr);
      costJb += 0.5 * ss[jj] * dot_product(*vvecs_[jj], *zvecs_[jj]) * ss[jj];
    }
    double costJoJc = costJ - costJb;

    // Gradient norm in precond metric --> sqrt(r'z) --> beta * s_{i}
    double rznorm = beta*std::abs(ss[jiter]);
    normReduction = rznorm/beta0;

    Log::info() << "DRPLanczos end of iteration " << jiter+1 << std::endl;
    printNormReduction(jiter+1, rznorm, normReduction);
    printQuadraticCostFunction(jiter+1, costJ, costJb, costJoJc);

    if (normReduction < tolerance) {
      Log::info() << "DRPLanczos: Achieved required reduction in residual norm." << std::endl;
      break;
    }
  }

  dx.zero();
  dxh.zero();

  // Calculate the solution (dxh = Binv dx)
  for (unsigned int jj = 0; jj < ss.size(); ++jj) {
    dx.axpy(ss[jj], *zvecs_[jj]);
    dxh.axpy(ss[jj], *hvecs_[jj]);
  }

  // Compute and save the eigenvectors
  writeEigenvectors(diagConf_, alphas_, betas_, dd, zvecs_, hvecs_, HtRinvH, pr, vv, zz);

  util::printRunStats("DRPLanczos end");
  return normReduction;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_DRPLANCZOSMINIMIZER_H_
