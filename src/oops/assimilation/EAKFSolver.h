/*
 * (C) Copyright 2026-2026 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_EAKFSOLVER_H_
#define OOPS_ASSIMILATION_EAKFSOLVER_H_

#include <algorithm>
#include <string>

#include "oops/assimilation/SequentialEnsembleSolver.h"

namespace oops {

/// This implements the Ensemble Adjustment Kalman Filter
template<typename MODEL, typename OBS>
class EAKFSolver : public SequentialEnsembleSolver<MODEL, OBS> {
  typedef Geometry<MODEL>   Geometry_;
  typedef ObsSpace<OBS>     ObsSpace_;
  typedef ObsSpaces<OBS>    ObsSpaces_;
  typedef StateSet<MODEL>   StateSet_;

 public:
  static const std::string classname() {return "oops::EAKFSolver";}

  EAKFSolver(ObsSpaces_ & obspaces, const Geometry_ & geometry,
             const eckit::Configuration & config, size_t nens,
             const StateSet_ & xbmean, const Variables & incvars)
    : SequentialEnsembleSolver<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars) {}

  virtual ~EAKFSolver() = default;

 protected:
  void obsEnsembleUpdate(const Eigen::VectorXd & yb_k,
                         const double omb_k,
                         const double oberr_variance_k,
                         Eigen::VectorXd & delta_y_k) override;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void EAKFSolver<MODEL, OBS>::obsEnsembleUpdate(const Eigen::VectorXd & yb_k,
                                               const double omb_k,
                                               const double oberr_variance_k,
                                               Eigen::VectorXd & delta_y_k) {
  // Implement a basic EAKF algorithm for updating the observation ensemble

  const Eigen::Index N = yb_k.size();
  const double yb_mean = yb_k.mean();

  // sample variance (unbiased)
  const double Pb = (yb_k.array() - yb_mean).square().sum() / (N - 1);
  if (Pb <= 0.0) {
    // hopefully we never get here, but just in case, exit gracefully
    delta_y_k.setZero();
    return;
  }

  // Kalman gain in observation space for scalar obs
  const double K = Pb / (Pb + oberr_variance_k);

  // posterior variance and perturbation scaling
  const double Pb_a = (1.0 - K) * Pb;
  const double scale = std::sqrt(std::max(0.0, Pb_a / Pb));

  // analysis mean shift
  const double mean_shift = K * omb_k;

  // fill in increments to the prior ensemble
  for (Eigen::Index i = 0; i < N; ++i) {
    const double pert = yb_k(i) - yb_mean;
    delta_y_k(i) = mean_shift + (scale - 1.0) * pert;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_EAKFSOLVER_H_
