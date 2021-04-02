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

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/gletkfInterface.h"
#include "oops/assimilation/LETKFSolver.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/ObsErrors.h"
#include "oops/util/Logger.h"

namespace oops {

/// Local Ensemble Tranform Kalman Filter solver using GSI GLETKF solver
template <typename MODEL, typename OBS>
class LETKFSolverGSI : public LETKFSolver<MODEL, OBS> {
  typedef Departures<OBS>           Departures_;
  typedef DeparturesEnsemble<OBS>   DeparturesEnsemble_;
  typedef Geometry<MODEL>           Geometry_;
  typedef ObsErrors<OBS>            ObsErrors_;
  typedef ObsSpaces<OBS>            ObsSpaces_;
 public:
  LETKFSolverGSI(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t);

  /// Computes weights
  void computeWeights(const Departures_ &, const DeparturesEnsemble_ &,
                      const Departures_ &) override;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
LETKFSolverGSI<MODEL, OBS>::LETKFSolverGSI(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                           const eckit::Configuration & config, size_t nens)
  : LETKFSolver<MODEL, OBS>(obspaces, geometry, config, nens)
{
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LETKFSolverGSI<MODEL, OBS>::computeWeights(const Departures_ & dy,
                                                const DeparturesEnsemble_ & Yb,
                                                const Departures_ & R_invvar) {
  // compute transformation matrix, save in Wa_, wa_
  // uses GSI GETKF code
  const int nobsl = dy.nobs();

  // cast oops objects to eigen<double> obejects
  // then cast eigen<doble> to eigen<float>
  Eigen::MatrixXd edy = dy.packEigen();
  Eigen::MatrixXf edy_f = edy.cast<float>();

  Eigen::MatrixXd eYb = Yb.packEigen();
  Eigen::MatrixXf eYb_f = eYb.cast<float>();

  Eigen::VectorXd eR = R_invvar.packEigen();
  Eigen::VectorXf eR_f = eR.cast<float>();

  Eigen::MatrixXf Wa_f(this->nens_, this->nens_);
  Eigen::VectorXf wa_f(this->nens_);

  // call into GSI interface to compute Wa and wa
  const int neigv = 1;
  const int getkf_inflation = 0;
  const int denkf = 0;
  const int getkf = 0;
  letkf_core_f90(nobsl, eYb_f.data(), eYb_f.data(), edy_f.data(),
                 wa_f.data(), Wa_f.data(),
                 eR_f.data(), this->nens_, neigv,
                 getkf_inflation, denkf, getkf);
  this->Wa_ = Wa_f.cast<double>();
  this->wa_ = wa_f.cast<double>();
}

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LETKFSOLVERGSI_H_
