/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_LETKFSOLVERGSI_H_
#define OOPS_ASSIMILATION_LETKFSOLVERGSI_H_

#include <Eigen/Dense>
#include <cfloat>
#include <vector>

#include "oops/assimilation/gletkfInterface.h"
#include "oops/assimilation/LETKFSolver.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/ObsErrors.h"

namespace oops {
  class Variables;

/// Local Ensemble Tranform Kalman Filter solver using GSI GLETKF solver
template <typename MODEL, typename OBS>
class LETKFSolverGSI : public LETKFSolver<MODEL, OBS> {
  typedef Departures<OBS>           Departures_;
  typedef DeparturesEnsemble<OBS>   DeparturesEnsemble_;
  typedef Geometry<MODEL>           Geometry_;
  typedef ObsErrors<OBS>            ObsErrors_;
  typedef ObsSpaces<OBS>            ObsSpaces_;
  typedef State4D<MODEL>            State4D_;
 public:
  LETKFSolverGSI(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t,
                 const State4D_ &, const Variables &);

  /// Computes weights for ensemble update with local observations
  /// \param[in] omb      Observation departures (nlocalobs)
  /// \param[in] Yb       Ensemble perturbations (nens, nlocalobs)
  /// \param[in] invvarR  Inverse of observation error variances (nlocalobs)
  virtual void computeWeights(const Eigen::VectorXd & omb, const Eigen::MatrixXd & Yb,
                              const Eigen::VectorXd & invvarR);
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
LETKFSolverGSI<MODEL, OBS>::LETKFSolverGSI(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                           const eckit::Configuration & config, size_t nens,
                                           const State4D_ & xbmean, const Variables & incvars)
  : LETKFSolver<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars)
{
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LETKFSolverGSI<MODEL, OBS>::computeWeights(const Eigen::VectorXd & dy,
                                                const Eigen::MatrixXd & Yb,
                                                const Eigen::VectorXd & R_invvar) {
  // compute transformation matrix, save in Wa_, wa_
  // uses GSI GETKF code
  const int nobsl = dy.size();
  const LocalEnsembleSolverInflationParameters & inflopt = this->options_.infl;
  const float infl = inflopt.mult;

  // cast eigen<double> to eigen<float>
  Eigen::VectorXf dy_f = dy.cast<float>();
  Eigen::MatrixXf Yb_f = Yb.cast<float>();
  Eigen::VectorXf R_invvar_f = R_invvar.cast<float>();

  Eigen::MatrixXf Wa_f(this->nens_, this->nens_);
  Eigen::VectorXf wa_f(this->nens_);

  // call into GSI interface to compute Wa and wa
  const int neigv = 1;
  const int getkf_inflation = 0;
  const int denkf = 0;
  const int getkf = 0;
  letkf_core_f90(nobsl, Yb_f.data(), Yb_f.data(), dy_f.data(),
                 wa_f.data(), Wa_f.data(),
                 R_invvar_f.data(), this->nens_, neigv,
                 getkf_inflation, denkf, getkf, infl);
  this->Wa_ = Wa_f.cast<double>();
  this->wa_ = wa_f.cast<double>();
}

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LETKFSOLVERGSI_H_
