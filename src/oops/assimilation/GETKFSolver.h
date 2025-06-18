/*
 * (C) Copyright 2020-2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_GETKFSOLVER_H_
#define OOPS_ASSIMILATION_GETKFSOLVER_H_

#include <Eigen/Dense>
#include <cfloat>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/ETKFLinearAlgebra.h"
#include "oops/assimilation/gletkfInterface.h"
#include "oops/assimilation/LocalEnsembleSolver.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/Geometry.h"
#include "oops/base/IncrementEnsemble4D.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/StateEnsemble4D.h"
#include "oops/base/StateSet.h"
#include "oops/generic/VerticalLocEV.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/Logger.h"
#include "oops/util/printRunStats.h"
#include "oops/util/Timer.h"

namespace oops {
  class Variables;

/*!
 * An implementation of the GETKF from Lei 2018 JAMES
 *
 * Lei, L., Whitaker, J. S., & Bishop, C. ( 2018). Improving assimilation
 * of radiance observations by implementing model space localization in an
 * ensemble Kalman filter. Journal of Advances in Modeling Earth Systems, 10,
 * 3221– 3232. https://doi.org/10.1029/2018MS001468
 */
template <typename MODEL, typename OBS>
class DeterministicGETKF : public LocalEnsembleSolver<MODEL, OBS> {
  typedef Departures<OBS>             Departures_;
  typedef DeparturesEnsemble<OBS>     DeparturesEnsemble_;
  typedef Geometry<MODEL>             Geometry_;
  typedef GeometryIterator<MODEL>     GeometryIterator_;
  typedef Increment<MODEL>            Increment_;
  typedef Increment4D<MODEL>          Increment4D_;
  typedef IncrementEnsemble4D<MODEL>  IncrementEnsemble4D_;
  typedef LinearModel<MODEL>          LinearModel_;
  typedef Model<MODEL>                Model_;
  typedef ModelAuxControl<MODEL>      ModelAux_;
  typedef ModelAuxIncrement<MODEL>    ModelAuxInc_;
  typedef ObsAuxControls<OBS>         ObsAux_;
  typedef ObsAuxIncrements<OBS>       ObsAuxInc_;
  typedef ObsDataVector<OBS, int>     ObsDataInt_;
  typedef ObsEnsemble<OBS>            ObsEnsemble_;
  typedef ObsErrors<OBS>              ObsErrors_;
  typedef Observations<OBS>           Observations_;
  typedef Observers<MODEL, OBS>       Observers_;
  typedef ObserversTLAD<MODEL, OBS>   ObserversTLAD_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef PseudoModelState4D<MODEL>   PseudoModel_;
  typedef PseudoLinearModelIncrement4D<MODEL> PseudoLinearModel_;
  typedef State<MODEL>                State_;
  typedef StateSet<MODEL>             StateSet_;
  typedef StateEnsemble4D<MODEL>      StateEnsemble4D_;
  typedef VerticalLocEV<MODEL>        VerticalLocEV_;

 public:
  static const std::string classname() {return "oops::DeterministicGETKF";}

  /// Constructor (allocates Wa, wa, HZb_,
  /// saves options from the config, computes VerticalLocEV_)
  DeterministicGETKF(ObsSpaces_ &,
                     const Geometry_ &,
                     const eckit::Configuration &,
                     size_t,
                     const StateSet_ &,
                     const Variables &);

  Observations_ computeHofX(const StateEnsemble4D_ &, size_t, bool) override;

  /// entire KF update (computeWeights+applyWeights) for a grid point GeometryIterator_
  void measurementUpdate(const Eigen::VectorXd &,
                         const Eigen::VectorXd &,
                         const Departures_ &,
                         const IncrementEnsemble4D_ &,
                         const GeometryIterator_ &,
                         IncrementEnsemble4D_ &) override;

 protected:
  Eigen::MatrixXd Wa_;  // transformation matrix for ens. perts. Xa_=Xf*Wa
  Eigen::VectorXd wa_;  // transformation matrix for ens. mean xa_=xf*wa
  size_t nens_;
  const Geometry_ & geometry_;
  VerticalLocEV_ vertloc_;
  size_t neig_;
  size_t nanal_;
  bool fortranETKF_;

  DeparturesEnsemble_ HZb_;

 private:
  /// Computes weights for ensemble update with local observations
  /// \param[in] omb      Observation departures (nlocalobs)
  /// \param[in] Yb       Ensemble perturbations for all the background memebers
  ///                     (nens*neig, nlocalobs)
  /// \param[in] YbOrig   Ensemble perturbations for the members to be updated (nens, nlocalobs)
  /// \param[in] invVarR  Inverse of observation error variances (nlocalobs)
  void computeWeights(const Eigen::VectorXd & omb,
                      const Eigen::MatrixXf & Yb,
                      const Eigen::MatrixXf & YbOrig,
                      const Eigen::VectorXd & invVarR);

  /// Applies weights and adds posterior inflation
  void applyWeights(const IncrementEnsemble4D_ &,
                    IncrementEnsemble4D_ &,
                    const GeometryIterator_ &);
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
DeterministicGETKF<MODEL, OBS>::DeterministicGETKF(ObsSpaces_ & obspaces,
                                                   const Geometry_ & geometry,
                                                   const eckit::Configuration & config,
                                                   size_t nens,
                                                   const StateSet_ & xbmean,
                                                   const Variables & incvars)
  : LocalEnsembleSolver<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars),
    nens_(nens), geometry_(geometry),
    vertloc_(config.getSubConfiguration("local ensemble DA.vertical localization"), xbmean[0],
    incvars), neig_(vertloc_.neig()), nanal_(neig_*nens_),
    fortranETKF_(config.getBool("local ensemble DA.fortran ETKF", true)),
    HZb_(obspaces, nanal_)
{
  // pre-allocate transformation matrices
  Wa_.resize(nanal_, nens);
  wa_.resize(nanal_);
}

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
Observations<OBS> DeterministicGETKF<MODEL, OBS>::computeHofX(const StateEnsemble4D_ & ens_xx,
                                                              size_t iteration, bool readFromFile) {
  util::Timer timer(classname(), "computeHofX");

  ModelAux_ moderr(geometry_, eckit::LocalConfiguration());
  ModelAuxInc_  moderrinc(geometry_, eckit::LocalConfiguration());
  ObsAux_  obsaux(this->obspaces_, this->observersconf_);
  ObsAuxInc_  obsauxinc(this->obspaces_, this->observersconf_);

  Observations_ yb_mean(this->obspaces_);

  if (readFromFile) {
    // compute/read H(x) for the original ensemble members
    // also computes omb_
    yb_mean = LocalEnsembleSolver<MODEL, OBS>::computeHofX(ens_xx, iteration, readFromFile);

    // read modulated ensemble
    Observations_ ytmp(yb_mean);
    size_t ii = 0;
    for (size_t iens = 0; iens < nens_; ++iens) {
      Log::info() << " DeterministicGETKF::computeHofX starting ensemble member "
                  << iens+1 << std::endl;
      util::printRunStats("DeterministicGETKF read hofx");
      for (size_t ieig = 0; ieig < neig_; ++ieig) {
        ytmp.read("hofxm"+std::to_string(iteration)+"_"+std::to_string(ieig+1)+
                      "_"+std::to_string(iens+1));
        HZb_.setData(ii, ytmp - yb_mean);
        ii = ii + 1;
      }
    }
  } else {
    const util::Duration default_tstep = (this->obspaces_.windowEnd()
                                        - this->obspaces_.windowStart()) * 2;
    const std::vector<util::DateTime> times = ens_xx[0].validTimes();
    const util::Duration flength = times[times.size()-1] - times[0];

    // Setup PseudoLinearModelIncrement4D to run on ensemble perturbation
    Increment4D_ dx(geometry_, this->incvars_, times);

    // modulate ensemble of obs
    IncrementEnsemble4D_ Ztmp(geometry_, this->incvars_, times, neig_);
    eckit::LocalConfiguration config;
    config.set("save hofx", false);
    config.set("save qc", false);
    config.set("save obs errors", false);
    config.set("iteration", std::to_string(iteration));
    size_t ii = 0;

    if (this->useLinearObserver()) {
      std::vector<ObsDataInt_> qcflags;
      for (size_t jj = 0; jj < this->obspaces_.size(); ++jj) {
        ObsDataInt_ qc(this->obspaces_[jj], this->obspaces_[jj].obsvariables());
        qcflags.push_back(qc);
      }
      // set up postprocessors for the linear model run on ensemble perturbations
      PostProcessor<Increment_> posttl;
      PostProcessorTLAD<MODEL> posttrajtl;

      this->computeHofX4D(config, this->xbmean_, yb_mean, flength, default_tstep, obsaux, moderr);

      // for linear H, yb_mean==y_mean_xb
      Observations_ y_mean_xb(yb_mean);

      // set QC for the mean
      config.set("save qc", true);
      config.set("save obs errors", true);
      y_mean_xb.save("hofx_y_mean_xb"+std::to_string(iteration));

      // QC flags and Obs errors are set to that of the H(mean(Xb))
      this->R_->save("ObsError");
      this->initializeAssimilatedMask();

      // mask H(x) ensemble perturbations - i.e. make sure that obs that have
      // failed QC on one ensemble member fail for all (this is for the case where
      // different QC procedures are done on different ensemble members)
      Departures_ tmpDeps(this->obspaces_);
      for (size_t iens = 0; iens < ens_xx.size(); ++iens) {
        tmpDeps.zero();
        tmpDeps = this->Yb_.getData(iens);
        this->updateAssimilatedMask(tmpDeps);
        this->applyAssimilatedMask(tmpDeps);
        this->Yb_.setData(iens, tmpDeps);
      }

      // calculate obs departures
      Observations_ yobs(this->obspaces_, "ObsValue");
      this->omb_ = yobs - yb_mean;
      // Need to mask out any missing departures as well as those that have failed QC
      this->updateAssimilatedMask(this->omb_);
      this->applyAssimilatedMask(this->omb_);

      // add linearized H(x) to the linear model postprocessor
      this->linear_hofx_->initializeTL(posttrajtl);

      for (size_t iens = 0; iens < ens_xx.size(); ++iens) {
        Log::info() << " DeterministicGETKF::computeHofX starting ensemble member "
                    << iens+1 << std::endl;
        util::printRunStats("DeterministicGETKF calculate hofx");
        Log::info() << " GETKFSolver::computeHofX starting ensemble member " << iens+1 << std::endl;
        util::printRunStats("GETKFSolver calculate hofx");
        tmpDeps.zero();
        // Setup PseudoLinearModelIncrement4D to run on ensemble perturbation
        dx.diff(ens_xx[iens], this->xbmean_);
        // Approximate H(x) using linearized model and linearized observer
        this->applyLinearToPerturbations(dx, flength, default_tstep, obsauxinc, moderrinc,
                                      posttl, posttrajtl, tmpDeps);
        (this->Yb_).setData(iens, tmpDeps);
        Observations_ tmpObs(yb_mean);
        tmpObs += this->Yb_.getData(iens);
        Log::test() << "H(x) for member " << iens+1 << ":" << std::endl << tmpObs << std::endl;
        tmpObs.save("hofx"+std::to_string(iteration)+"_"+std::to_string(iens+1));

        // observe modulated members
        vertloc_.modulateIncrement(dx, Ztmp);
        for (size_t ieig = 0; ieig < neig_; ++ieig) {
          this->applyLinearToPerturbations(Ztmp[ieig], flength, default_tstep, obsauxinc, moderrinc,
                                           posttl, posttrajtl, tmpDeps);
          HZb_.setData(ii, tmpDeps);
          Observations_ tmpObs(yb_mean);
          tmpObs += HZb_.getData(ii);
          tmpObs.save("hofxm"+std::to_string(iteration)+"_"+std::to_string(ieig+1)+
                        "_"+std::to_string(iens+1));
          ii = ii + 1;
        }
      }
    } else {
      // compute/read H(x) for the original ensemble members
      // also computes omb_
      yb_mean = LocalEnsembleSolver<MODEL, OBS>::computeHofX(ens_xx, iteration, readFromFile);

      for (size_t iens = 0; iens < nens_; ++iens) {
        Log::info() << " DeterministicGETKF::computeHofX starting ensemble member "
                    << iens+1 << std::endl;
        util::printRunStats("DeterministicGETKF calculate hofx");
        dx.diff(ens_xx[iens], this->xbmean_);
        vertloc_.modulateIncrement(dx, Ztmp);
        for (size_t ieig = 0; ieig < neig_; ++ieig) {
          StateSet_ tmpState = this->xbmean_;
          tmpState += Ztmp[ieig];
          Observations_ tmpObs(this->obspaces_);
          this->computeHofX4D(config, tmpState, tmpObs, flength, default_tstep, obsaux, moderr);
          HZb_.setData(ii, tmpObs - yb_mean);
          tmpObs.save("hofxm"+std::to_string(iteration)+"_"+std::to_string(ieig+1)+
                        "_"+std::to_string(iens+1));
          ii = ii + 1;
        }
      }
    }
  }
  // Update mask again, this time for the modulated ensemble members
  Departures_ tmpDeps(this->obspaces_);
  for (size_t iens = 0; iens < nanal_; ++iens) {
    tmpDeps.zero();
    tmpDeps = this->HZb_.getData(iens);
    this->updateAssimilatedMask(tmpDeps);
    this->applyAssimilatedMask(tmpDeps);
    this->HZb_.setData(iens, tmpDeps);
  }
  return yb_mean;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void DeterministicGETKF<MODEL, OBS>::computeWeights(const Eigen::VectorXd & dy,
                                                    const Eigen::MatrixXf & Yb,
                                                    const Eigen::MatrixXf & YbOrig,
                                                    const Eigen::VectorXd & invVarR) {
  // compute transformation matrix, save in Wa_, wa_
  // Yb(nobs,neig*nens), YbOrig(nobs,nens)
  util::Timer timer(classname(), "computeWeights");
  const LocalEnsembleSolverInflationParameters & inflopt = this->options_.infl;
  const float infl = inflopt.mult;

  Eigen::MatrixXf Wa_f(nanal_, this->nens_);
  Eigen::VectorXf wa_f(nanal_);

  if (fortranETKF_) {
    Eigen::MatrixXf Wa_f(nanal_, this->nens_);
    Eigen::VectorXf wa_f(nanal_);

    // cast eigen<double> to eigen<float>
    const Eigen::VectorXf dy_f = dy.cast<float>();
    const Eigen::VectorXf invVarR_f = invVarR.cast<float>();

    // call into GSI interface to compute Wa and wa
    const int nobsl = dy.size();
    const int getkf_inflation = 0;
    const int denkf = 0;
    const int getkf = 1;
    letkf_core_f90(nobsl, Yb.data(), YbOrig.data(), dy_f.data(),
                   wa_f.data(), Wa_f.data(),
                   invVarR_f.data(), nanal_, neig_,
                   getkf_inflation, denkf, getkf, infl);

    this->Wa_ = Wa_f.cast<double>();
    this->wa_ = wa_f.cast<double>();
  } else {
    oops::detGETKF_computeWeights(dy, Yb, YbOrig, invVarR, infl, wa_, Wa_);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void DeterministicGETKF<MODEL, OBS>::applyWeights(const IncrementEnsemble4D_ & bkg_pert,
                                                  IncrementEnsemble4D_ & ana_pert,
                                                  const GeometryIterator_ & i) {
  util::Timer timer(classname(), "applyWeights");

  oops::detETKF_applyWeights<MODEL>(bkg_pert,
                                    ana_pert,
                                    i,
                                    this->wa_,
                                    this->Wa_,
                                    this->options_.infl,
                                    &(vertloc_));
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void DeterministicGETKF<MODEL, OBS>::measurementUpdate(const Eigen::VectorXd & local_omb_vec,
                                                       const Eigen::VectorXd & local_invVarR_vec,
                                                       const Departures_ & locvector,
                                                       const IncrementEnsemble4D_ & bkg_pert,
                                                       const GeometryIterator_ & i,
                                                       IncrementEnsemble4D_ & ana_pert) {
  const Eigen::MatrixXf local_Yb_mat_f = this->Yb_.packEigen(locvector);
  const Eigen::MatrixXf local_HZ_mat_f = this->HZb_.packEigen(locvector);
  this->computeWeights(local_omb_vec, local_HZ_mat_f, local_Yb_mat_f, local_invVarR_vec);
  this->applyWeights(bkg_pert, ana_pert, i);
}

}  // namespace oops
#endif  // OOPS_ASSIMILATION_GETKFSOLVER_H_
