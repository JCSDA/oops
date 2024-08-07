/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_QUADRATURESOLVER_H_
#define OOPS_ASSIMILATION_QUADRATURESOLVER_H_

#include <cfloat>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/QuadratureRules.h"
#include "oops/assimilation/HMatrix.h"
#include "oops/assimilation/HtMatrix.h"
#include "oops/assimilation/BMatrix.h"
#include "oops/assimilation/RinvSqrtMatrix.h"
#include "oops/assimilation/NormalizedHBHtMatrix.h"
#include "oops/assimilation/SLCG.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL, typename OBS>
class QuadratureSolver {
  typedef CostFunction<MODEL, OBS>          CostFct_;
  typedef ControlIncrement<MODEL, OBS>      CtrlInc_;
  typedef DualVector<MODEL, OBS>            Dual_;
  typedef HMatrix<MODEL, OBS>               H_;
  typedef HtMatrix<MODEL, OBS>              Ht_;
  typedef BMatrix<MODEL, OBS>               B_;
  typedef RinvSqrtMatrix<MODEL, OBS>        R_invsqrt_;
  typedef NormalizedHBHtMatrix<MODEL, OBS>  NormalHBHt_;

  public:
    QuadratureSolver(const CostFct_ & J);
    ~QuadratureSolver() {}
    void solve(CtrlInc_ dx, const eckit::Configuration & config);

  private:
    const H_          H_mat_;
    const Ht_         Ht_mat_;
    const B_          B_mat_;
    const R_invsqrt_  R_invsqrt_mat_;
    const NormalHBHt_ NormalHBHt_mat_;
    const double      PI;
};

// =============================================================================

template<typename MODEL, typename OBS>
QuadratureSolver<MODEL, OBS>::QuadratureSolver(const CostFct_ & J)
  : H_mat_(J), Ht_mat_(J), B_mat_(J), R_invsqrt_mat_(J), NormalHBHt_mat_(J), PI(4*atan(1)) {}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void QuadratureSolver<MODEL, OBS>::solve(CtrlInc_ dx, const eckit::Configuration & config) {
  Log::info() << "QuadratureSolver: Starting." << std::endl;

  int quadsize     = config.getInt("quadsize");
  int maxiters     = config.getInt("maxiters");
  double tolerance = config.getDouble("tolerance");

  Log::info() << "QuadratureSolver:" << std::endl;
  Log::info() << "  quadrature size           = " << quadsize << std::endl;
  Log::info() << "  max iteration count       = " << maxiters << std::endl;
  Log::info() << "  linear solution tolerance = " << tolerance << std::endl;
  Log::info() << "QuadratureSolver: Mapping increment to observation space." << std::endl;

  Dual_ dy, dz_in;
  H_mat_.multiply(dx, dy);
  R_invsqrt_mat_.multiply(dy, dz_in);

// Set up nodes and weights for Gauss-Legendre quadrature
  std::vector<double> weights;
  std::vector<double> nodes;
  gaussLegendre(quadsize, nodes, weights);
  prepareEAKFQuad(quadsize, nodes, weights);

  Log::info() << "QuadratureSolver: Beginning linear system solves." << std::endl;
  std::vector<Dual_> dz_sols;
  SLCG(dz_sols, NormalHBHt_mat_, dz_in, nodes, quadsize, maxiters, tolerance);

  Log::info() << "QuadratureSolver: Recombining solutions." << std::endl;

  Dual_ dz_out(dz_sols[0]);
  dz_out.zero();

  for(int q = 0; q < quadsize; q++) {
    dz_out.axpy(weights[q], dz_sols[q]);
  }

  Log::info() << "QuadratureSolver: Mapping increment to state space." << std::endl;

  Dual_    tmp_dual(dz_out);
  CtrlInc_ tmp_ctrl(dx), dx_inc(dx);

  R_invsqrt_mat_.multiply(dz_out, tmp_dual);
  Ht_mat_.multiply(tmp_dual, tmp_ctrl);
  B_mat_.multiply(tmp_ctrl, dx_inc);

  Log::info() << "QuadratureSolver: Applying adjustment to input increment." << std::endl;
  dx.axpy(-1, dx_inc);

  Log::info() << "QuadratureSolver: Finished." << std::endl;
}

}  // namespace oops
#endif  // OOPS_ASSIMILATION_QUADRATURESOLVER_H_
