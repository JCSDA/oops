/*
 * (C) Copyright 2020-2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#ifndef OOPS_ASSIMILATION_LETKFSOLVER_H_
#define OOPS_ASSIMILATION_LETKFSOLVER_H_

#include <Eigen/Dense>
#include <cfloat>
#include <memory>
#include <string>
#include <vector>

#include "oops/assimilation/ETKFLinearAlgebra.h"
#include "oops/assimilation/LocalEnsembleSolver.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/Geometry.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/ObsLocalizations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/StateSet.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/util/Logger.h"

namespace oops {
  class Variables;

/// Local Ensemble Tranform Kalman Filter solver
/*!
 * An implementation of the LETKF from Hunt et al. 2007
 * this version is implemented using Eigen algebra and
 * temporary Eigen matrices for Xa and Xb
 * this verion implements RTPP and RTPS.
 *
 * Hunt, B. R., Kostelich, E. J., & Szunyogh, I. (2007). Efficient data
 * assimilation for spatiotemporal chaos: A local ensemble transform Kalman
 * filter. Physica D: Nonlinear Phenomena, 230(1-2), 112-126.
 */
template <typename MODEL, typename OBS>
class DeterministicLETKF : public LocalEnsembleSolver<MODEL, OBS> {
  typedef Departures<OBS>             Departures_;
  typedef DeparturesEnsemble<OBS>     DeparturesEnsemble_;
  typedef Geometry<MODEL>             Geometry_;
  typedef GeometryIterator<MODEL>     GeometryIterator_;
  typedef IncrementSet<MODEL>         IncrementSet_;
  typedef ObsErrors<OBS>              ObsErrors_;
  typedef ObsLocalizations<MODEL, OBS> ObsLocalizations_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef StateSet<MODEL>             StateSet_;

 public:
  static const std::string classname() {return "oops::DeterministicLETKF";}

  DeterministicLETKF(ObsSpaces_ &,
                     const Geometry_ &,
                     const eckit::Configuration &,
                     size_t,
                     const StateSet_ &,
                     const Variables &);

  /// KF update + posterior inflation at a grid point location (GeometryIterator_)
  void measurementUpdate(const Eigen::VectorXd &,
                         const Eigen::VectorXd &,
                         const Departures_ &,
                         const IncrementSet_ &,
                         const GeometryIterator_ &,
                         IncrementSet_ &) override;

 protected:
  /// Computes weights for ensemble update with local observations
  /// \param[in] omb      Observation departures (nlocalobs)
  /// \param[in] Yb       Ensemble perturbations (nens, nlocalobs)
  /// \param[in] invVarR  Inverse of observation error variances (nlocalobs)

  void computeWeights(const Eigen::VectorXd & omb,
                              const Eigen::MatrixXf & Yb,
                              const Eigen::VectorXd & invVarR);

  /// Applies weights and adds posterior inflation
  virtual void applyWeights(const IncrementSet_ &,
                            IncrementSet_ &,
                            const GeometryIterator_ &);

  Eigen::MatrixXd Wa_;  // transformation matrix for ens. perts. Xa=Xf*Wa
  Eigen::VectorXd wa_;  // transformation matrix for ens. mean xa=xf*wa

  const size_t nens_;   // ensemble size
  bool fortranETKF_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
DeterministicLETKF<MODEL, OBS>::DeterministicLETKF(ObsSpaces_ & obspaces,
                                                   const Geometry_ & geometry,
                                                   const eckit::Configuration & config,
                                                   size_t nens,
                                                   const StateSet_ & xbmean,
                                                   const Variables & incvars)
  : LocalEnsembleSolver<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars),
    nens_(nens),
    fortranETKF_(config.getBool("local ensemble DA.fortran ETKF", false))
{
  Log::trace() << "DeterministicLETKF<MODEL, OBS>::create starting" << std::endl;
  Log::info() << "Using EIGEN implementation of LETKF" << std::endl;

  // pre-allocate transformation matrices
  Wa_.resize(nens_, nens_);
  wa_.resize(nens_);

  Log::trace() << "DeterministicLETKF<MODEL, OBS>::create done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void DeterministicLETKF<MODEL, OBS>::measurementUpdate(const Eigen::VectorXd & local_omb_vec,
                                                       const Eigen::VectorXd & local_invVarR_vec,
                                                       const Departures_ & locvector,
                                                       const IncrementSet_ & bkg_pert,
                                                       const GeometryIterator_ & i,
                                                       IncrementSet_ & ana_pert) {
  const Eigen::MatrixXf local_Yb_mat_f = this->Yb_.packEigen(locvector);
  this->computeWeights(local_omb_vec, local_Yb_mat_f, local_invVarR_vec);
  this->applyWeights(bkg_pert, ana_pert, i);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void DeterministicLETKF<MODEL, OBS>::computeWeights(const Eigen::VectorXd & dy,
                                                    const Eigen::MatrixXf & Yb,
                                                    const Eigen::VectorXd & invVarR ) {
  // compute transformation matrix, save in Wa_, wa_
  // uses C++ eigen interface
  // implements LETKF from Hunt et al. 2007
  util::Timer timer(classname(), "computeWeights");

  const LocalEnsembleSolverInflationParameters & inflopt = this->options_.infl;
  const double infl = inflopt.mult;

  if (fortranETKF_) {
    Eigen::MatrixXf Wa_f(this->nens_, this->nens_);
    Eigen::VectorXf wa_f(this->nens_);

    // cast eigen<double> to eigen<float>
    const Eigen::VectorXf dy_f = dy.cast<float>();
    const Eigen::VectorXf invVarR_f = invVarR.cast<float>();

    // call into GSI interface to compute Wa and wa
    const int nobsl = dy.size();
    const int getkf_inflation = 0;
    const int denkf = 0;
    const int getkf = 0;
    const int neigv = 1;
    letkf_core_f90(nobsl, Yb.data(), Yb.data(), dy_f.data(),
                   wa_f.data(), Wa_f.data(),
                   invVarR_f.data(), this->nens_, neigv,
                   getkf_inflation, denkf, getkf, infl);

    this->Wa_ = Wa_f.cast<double>();
    this->wa_ = wa_f.cast<double>();
  } else {
    oops::detLETKF_computeWeights(dy, Yb, invVarR, (nens_ - 1) / infl, wa_, Wa_);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void DeterministicLETKF<MODEL, OBS>::applyWeights(const IncrementSet_ & bkg_pert,
                                                  IncrementSet_ & ana_pert,
                                                  const GeometryIterator_ & i) {
  util::Timer timer(classname(), "applyWeights");

  oops::detETKF_applyWeights<MODEL>(bkg_pert,
                                    ana_pert,
                                    i,
                                    this->wa_,
                                    this->Wa_,
                                    this->options_.infl);
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LETKFSOLVER_H_
