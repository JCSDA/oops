/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_RPCGMINIMIZER_H_
#define OOPS_ASSIMILATION_RPCGMINIMIZER_H_

#include <string>
#include <vector>

#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualMinimizer.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/HBHtMatrix.h"
#include "oops/assimilation/RinvMatrix.h"
#include "oops/base/IdentityMatrix.h"
#include "util/dot_product.h"
#include "util/formats.h"
#include "util/Logger.h"

namespace oops {

/// RPCG Minimizer
/*!
 * \brief Augmented Restricted Preconditioned Conjugate Gradients.
 *
 * This solver is based on the algorithm proposed in Gratton and Tshimanga,
 * QJRMS, 135: 1573-1585 (2009). It performs minimization in observation space.
 *
 * It solves the linear system \f$ (I + R^{-1}HBH^T) lambda  = H^T R^{-1}d \f$
 * with \f$ HBH^T \f$ inner-product in the augmented observation space.
 *
 * Note that the traditional \f$ B\f$-preconditioning in model space
 * corresponds to \f$I\f$ for this algorithm.
 *
 * A second-level preconditioner, \f$ G \f$, must be symmetric and
 * positive definite with respect to \f$ HBH^T \f$ inner-product.
 * Possible preconditioning is detailed in S. Gurol, PhD Manuscript, 2013.
 *
 * On entry:
 * -    vv      =  starting point, \f$ v_0 \f$.
 * -    vvp     =  starting point augmented part, \f$ vp_0 \f$
 * -    rr      = right hand side.
 * -    dy      = \f$ H (xb - xk) \f$
 * -    sigma   = \f$ (xb - xk)^T B^{-1} (xb - xk)\f$.
 * -    HBHt    = \f$ HBH^T \f$.
 * -    Rinv    = \f$ R^{-1} \f$.
 * -    G       = preconditioner \f$ G \f$.
 * -    Gt      = preconditioner transpose \f$ G^T \f$.
 *
 * On exit, vv and vvp will contain the solution \f$ lambda = [vv; vp] \f$
 *
 *  The return value is the achieved reduction in preconditioned residual norm.
 *
 *  Iteration will stop if the maximum iteration limit "maxiter" is reached
 *  or if the preconditioned residual norm reduces by a factor of "tolerance".
 */

// -----------------------------------------------------------------------------

template<typename MODEL> class RPCGMinimizer : public DualMinimizer<MODEL> {
  typedef CostFunction<MODEL>        CostFct_;
  typedef DualVector<MODEL>          Dual_;
  typedef HBHtMatrix<MODEL>          HBHt_;
  typedef RinvMatrix<MODEL>          Rinv_;

 public:
  const std::string classname() const override {return "RPCGMinimizer";}
  RPCGMinimizer(const eckit::Configuration &, const CostFct_ & J): DualMinimizer<MODEL>(J) {}
  ~RPCGMinimizer() {}

 private:
  double solve(Dual_ &, double &, Dual_ &, const HBHt_ &, const Rinv_ &,
               const int &, const double &, Dual_ &, const double &) override;
};

// =============================================================================

template<typename MODEL>
double RPCGMinimizer<MODEL>::solve(Dual_ & vv, double & vvp, Dual_ & rr,
                                   const HBHt_ & HBHt, const Rinv_ & Rinv,
                                   const int & maxiter, const double & tolerance,
                                   Dual_ & dy, const double & sigma) {
  IdentityMatrix<Dual_> precond;
  IdentityMatrix<Dual_> precondt;

  Dual_ zz(vv);
  Dual_ ww(vv);
  Dual_ pp(vv);
  Dual_ tt(vv);
  Dual_ ll(vv);
  Dual_ qq(vv);
  Dual_ w(vv);
  Dual_ v(vv);

  // augmented part of the vectors
  double rrp = 1.0;
  double llp;
  double wwp;
  double zzp;
  double ppp = 0.0;
  double ttp = 0.0;
  double qqp;
  double vp = 1.0;
  double wp;

  std::vector<Dual_> vVEC;  // required for re-orthogonalization
  std::vector<Dual_> wVEC;  // required for re-orthogonalization
  std::vector<double> vpVEC;
  std::vector<double> wpVEC;

  // llaug = HBHtaug rraug
  HBHt.multiply(rr, ll);
  ll.axpy(rrp, dy);
  llp = dot_product(dy, rr) + sigma*rrp;

  // wwaug = Gtaug llaug
  precondt.multiply(ll, ww);  // UPDATE WHEN G IS NOT IDENTITY
  wwp = llp;

  // zzaug = Gaug  rraug
  precond.multiply(rr, zz);   // UPDATE WHEN G IS NOT IDENTITY
  zzp = rrp;

  double dotwr = dot_product(rr, ww) + rrp*wwp;
  double normReduction = 1.0;
  double dotw0r0 = dotwr;
  double dotwr_old = 0.0;

  v = rr;
  v *= 1/sqrt(dotwr);
  vp = rrp;
  vp *= 1/sqrt(dotwr);
  w = ww;
  w *= 1/sqrt(dotwr);
  wp = wwp;
  wp *= 1/sqrt(dotwr);

  vVEC.clear();
  vpVEC.clear();
  wVEC.clear();
  wpVEC.clear();

  vVEC.push_back(v);
  vpVEC.push_back(vp);
  wVEC.push_back(w);
  wpVEC.push_back(wp);

  Log::info() << std::endl;
  for (int jiter = 0; jiter < maxiter; ++jiter) {
    Log::info() << "RPCG Starting Iteration " << jiter+1 << std::endl;

    if (jiter > 0) {
      double beta = dotwr/dotwr_old;
      Log::debug() << "RPCG beta = " << beta << std::endl;

      pp *= beta;
      ppp *= beta;

      tt *= beta;
      ttp *= beta;
    }

    pp += zz;
    ppp += zzp;  // ppaug = zzaug + beta*ppaug

    tt += ww;
    ttp += wwp;  // ttaug = wwaug + beta*ttaug

    // (RinvHBHt + I) pp
    Rinv.multiply(tt, qq);
    qq += pp;
    qqp = ppp;  // qqaug = ppaug + Rinv_aug ttaug

    double alpha = dotwr/(dot_product(qq, tt) + qqp * ttp);
    Log::debug() << "RPCG alpha = " << alpha << std::endl;

    vv.axpy(alpha, pp);
    vvp = vvp + alpha*ppp;  // vvaug = vvaug + alpha*ppaug

    rr.axpy(-alpha, qq);
    rrp = rrp - alpha*qqp;  // rraug = rraug - alpha*qqaug

    // Re-orthogonalization
    for (int iiter = 0; iiter < jiter; ++iiter) {
      double proj = dot_product(rr, wVEC[iiter]) + rrp * wpVEC[iiter];
      rr.axpy(-proj, vVEC[iiter]);
      rrp -= proj*vpVEC[iiter];
    }

    // llaug = HBHt_aug rraug
    HBHt.multiply(rr, ll);
    ll.axpy(rrp, dy);
    llp = dot_product(dy, rr) + sigma * rrp;

    // wwaug = Gt_aug llaug
    precondt.multiply(ll, ww);  // UPDATE WHEN G IS NOT IDENTITY
    wwp = llp;

    // zzaug = Gaug  rraug
    precond.multiply(rr, zz);   // UPDATE WHEN G IS NOT IDENTITY
    zzp = rrp;

    dotwr_old = dotwr;
    dotwr = dot_product(ww, rr) + rrp * wwp;

    v = rr;
    v *= 1/sqrt(dotwr);
    vp = rrp;
    vp *= 1/sqrt(dotwr);
    w = ww;
    w *= 1/sqrt(dotwr);
    wp = wwp;
    wp *= 1/sqrt(dotwr);

    vVEC.push_back(v);
    vpVEC.push_back(vp);
    wVEC.push_back(w);
    wpVEC.push_back(wp);

    normReduction = sqrt(dotwr/dotw0r0);

    Log::info() << "RPCG end of iteration " << jiter+1 << ". Norm reduction= "
                << util::full_precision(normReduction) << std::endl << std::endl;

    if (normReduction < tolerance) {
      Log::info() << "RPCG: Achieved required reduction in residual norm." << std::endl;
      break;
    }
  }
  return normReduction;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_RPCGMINIMIZER_H_
