/*
 * (C) Copyright 2020 UCAR.
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
#include "oops/base/State4D.h"
#include "oops/base/StateEnsemble4D.h"
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
class GETKFSolver : public LocalEnsembleSolver<MODEL, OBS> {
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
  typedef State4D<MODEL>              State4D_;
  typedef StateEnsemble4D<MODEL>      StateEnsemble4D_;
  typedef VerticalLocEV<MODEL>        VerticalLocEV_;

 public:
  static const std::string classname() {return "oops::GETKFSolver";}

  /// Constructor (allocates Wa, wa, HZb_,
  /// saves options from the config, computes VerticalLocEV_)
  GETKFSolver(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t,
              const State4D_ &, const Variables &);

  Observations_ computeHofX(const StateEnsemble4D_ &, size_t, bool) override;

  /// entire KF update (computeWeights+applyWeights) for a grid point GeometryIterator_
  void measurementUpdate(const IncrementEnsemble4D_ &, const GeometryIterator_ &,
                         IncrementEnsemble4D_ &) override;

 private:
  /// Computes weights for ensemble update with local observations
  /// \param[in] omb      Observation departures (nlocalobs)
  /// \param[in] Yb       Ensemble perturbations for all the background memebers
  ///                     (nens*neig, nlocalobs)
  /// \param[in] YbOrig   Ensemble perturbations for the members to be updated (nens, nlocalobs)
  /// \param[in] invvarR  Inverse of observation error variances (nlocalobs)
  void computeWeights(const Eigen::VectorXd & omb, const Eigen::MatrixXd & Yb,
                      const Eigen::MatrixXd & YbOrig, const Eigen::VectorXd & invvarR);

  /// Applies weights and adds posterior inflation
  void applyWeights(const IncrementEnsemble4D_ &, IncrementEnsemble4D_ &,
                    const GeometryIterator_ &);

 private:
  // parameters
  size_t nens_;
  const Geometry_ & geometry_;
  VerticalLocEV_ vertloc_;
  size_t neig_;
  size_t nanal_;

  DeparturesEnsemble_ HZb_;

  Eigen::MatrixXd Wa_;  // transformation matrix for ens. perts. Xa_=Xf*Wa
  Eigen::VectorXd wa_;  // transformation matrix for ens. mean xa_=xf*wa
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GETKFSolver<MODEL, OBS>::GETKFSolver(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                const eckit::Configuration & config, size_t nens,
                                const State4D_ & xbmean, const Variables & incvars)
  : LocalEnsembleSolver<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars),
    nens_(nens), geometry_(geometry),
    vertloc_(config.getSubConfiguration("local ensemble DA.vertical localization"), xbmean[0],
    incvars), neig_(vertloc_.neig()), nanal_(neig_*nens_), HZb_(obspaces, nanal_)
{
  // pre-allocate transformation matrices
  Wa_.resize(nanal_, nens);
  wa_.resize(nanal_);
}

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
Observations<OBS> GETKFSolver<MODEL, OBS>::computeHofX(const StateEnsemble4D_ & ens_xx,
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
      Log::info() << " GETKFSolver::computeHofX starting ensemble member " << iens+1 << std::endl;
      util::printRunStats("GETKFSolver read hofx");
      for (size_t ieig = 0; ieig < neig_; ++ieig) {
        ytmp.read("hofxm"+std::to_string(iteration)+"_"+std::to_string(ieig+1)+
                      "_"+std::to_string(iens+1));
        HZb_[ii] = ytmp - yb_mean;
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
      this->R_.reset(new ObsErrors_(this->observersconf_, this->obspaces_));

      // Setup pseudo model to run on ensemble mean
      State_ init_xx = this->xbmean_[0];
      std::unique_ptr<PseudoModel_> pseudomodel(new PseudoModel_(this->xbmean_, default_tstep));
      const Model_ model(std::move(pseudomodel));

      // setup postprocessors and nonlinear observers for the "nonlinear" model run on the mean
      PostProcessor<State_> post;
      PostProcessorTLAD<MODEL> posttraj;
      Observers_ hofx(this->obspaces_, this->obsconf_);

      // setup postprocessors and linear observers for the "linear" model run on the ensemble
      // perturbations
      PostProcessor<Increment_> posttl;
      PostProcessorTLAD<MODEL> posttrajtl;
      ObserversTLAD_ linear_hofx(this->obspaces_, this->obsconf_);

      // initialize nonlinear model postprocessor
      hofx.initialize(this->geometry_, obsaux, *this->R_, post, config);

      // add linearized H(x) to the nonlinear model postprocessor
      linear_hofx.initializeTraj(this->geometry_, obsaux, posttraj);
      // create TrajectorySaver with hofx_linear, and enroll in post
      post.enrollProcessor(new TrajectorySaver<MODEL>(eckit::LocalConfiguration(),
                                                  this->geometry_, posttraj));
      // run nonlinear model on the ensemble mean
      model.forecast(init_xx, moderr, flength, post);
      // compute nonlinear H(x_mean)
      std::vector<ObsDataInt_> qcflags;
      for (size_t jj = 0; jj < this->obspaces_.size(); ++jj) {
        ObsDataInt_ qc(this->obspaces_[jj], this->obspaces_[jj].obsvariables());
        qcflags.push_back(qc);
      }
      hofx.finalize(yb_mean, qcflags);
      linear_hofx.finalizeTraj(qcflags);

      // for linear H, yb_mean==y_mean_xb
      Observations_ y_mean_xb(yb_mean);

      // set QC for the mean
      config.set("save qc", true);
      config.set("save obs errors", true);
      y_mean_xb.save("hofx_y_mean_xb"+std::to_string(iteration));

      // QC flags and Obs errors are set to that of the H(mean(Xb))
      this->R_->save("ObsError");
      // set inverse variances
      this->invVarR_.reset(new Departures_(this->R_->inverseVariance()));

      // mask H(x) ensemble perturbations
      for (size_t iens = 0; iens < ens_xx.size(); ++iens) {
        this->invVarR_->mask(this->Yb_[iens]);
        this->Yb_[iens].mask(*this->invVarR_);
      }

      // calculate obs departures and mask with qc flag
      Observations_ yobs(this->obspaces_, "ObsValue");
      this->omb_ = yobs - yb_mean;
      this->invVarR_->mask(this->omb_);
      this->omb_.mask(*this->invVarR_);

      // add linearized H(x) to the linear model postprocessor
      linear_hofx.initializeTL(posttrajtl);
      for (size_t iens = 0; iens < ens_xx.size(); ++iens) {
        Log::info() << " GETKFSolver::computeHofX starting ensemble member " << iens+1 << std::endl;
        util::printRunStats("GETKFSolver calculate hofx");

        dx.diff(ens_xx[iens], this->xbmean_);

        // observe original member
        Increment_ init_dx = dx[0];
        std::unique_ptr<PseudoLinearModel_> pseudolinearmodel =
           std::make_unique<PseudoLinearModel_>(dx, default_tstep);
        const LinearModel_ linear_model(std::move(pseudolinearmodel));
        // run linear model on the ensemble perturbation, compute linear H*dx
        linear_model.forecastTL(init_dx, moderrinc, flength, posttl, posttrajtl);
        linear_hofx.finalizeTL(obsauxinc, this->Yb_[iens]);

        Observations_ tmpObs(yb_mean);
        tmpObs += this->Yb_[iens];
        Log::test() << "H(x) for member " << iens+1 << ":" << std::endl << tmpObs << std::endl;
        tmpObs.save("hofx"+std::to_string(iteration)+"_"+std::to_string(iens+1));

        // observe modulated members
        vertloc_.modulateIncrement(dx, Ztmp);
        for (size_t ieig = 0; ieig < neig_; ++ieig) {
          std::unique_ptr<PseudoLinearModel_> pseudolinearmodel =
               std::make_unique<PseudoLinearModel_>(Ztmp[ieig], default_tstep);
          const LinearModel_ linear_model(std::move(pseudolinearmodel));
          // run linear model on the ensemble perturbation, compute linear H*dx
          Increment_ init_dx = Ztmp[ieig][0];
          linear_model.forecastTL(init_dx, moderrinc, flength, posttl, posttrajtl);
          linear_hofx.finalizeTL(obsauxinc, HZb_[ii]);
          Observations_ tmpObs(yb_mean);
          tmpObs += HZb_[ii];
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
        Log::info() << " GETKFSolver::computeHofX starting ensemble member " << iens+1 << std::endl;
        util::printRunStats("GETKFSolver calculate hofx");
        dx.diff(ens_xx[iens], this->xbmean_);
        vertloc_.modulateIncrement(dx, Ztmp);
        for (size_t ieig = 0; ieig < neig_; ++ieig) {
          State4D_ tmpState = this->xbmean_;
          tmpState += Ztmp[ieig];
          Observations_ tmpObs(this->obspaces_);
          this->computeHofX4DNonLinear(config, tmpState, tmpObs);
          HZb_[ii] = tmpObs - yb_mean;
          tmpObs.save("hofxm"+std::to_string(iteration)+"_"+std::to_string(ieig+1)+
                        "_"+std::to_string(iens+1));
          ii = ii + 1;
        }
      }
    }
  }
  return yb_mean;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GETKFSolver<MODEL, OBS>::computeWeights(const Eigen::VectorXd & dy,
                                             const Eigen::MatrixXd & Yb,
                                             const Eigen::MatrixXd & YbOrig,
                                             const Eigen::VectorXd & R_invvar) {
  // compute transformation matrix, save in Wa_, wa_
  // Yb(nobs,neig*nens), YbOrig(nobs,nens)
  // uses GSI GETKF code
  util::Timer timer(classname(), "computeWeights");
  const LocalEnsembleSolverInflationParameters & inflopt = this->options_.infl;
  const float infl = inflopt.mult;

  const int nobsl = dy.size();

  // cast eigen<double> to eigen<float>
  Eigen::VectorXf dy_f = dy.cast<float>();
  Eigen::MatrixXf Yb_f = Yb.cast<float>();
  Eigen::MatrixXf YbOrig_f = YbOrig.cast<float>();
  Eigen::VectorXf R_invvar_f = R_invvar.cast<float>();

  Eigen::MatrixXf Wa_f(nanal_, this->nens_);
  Eigen::VectorXf wa_f(nanal_);

  // call into GSI interface to compute Wa and wa
  const int getkf_inflation = 0;
  const int denkf = 0;
  const int getkf = 1;
  letkf_core_f90(nobsl, Yb_f.data(), YbOrig_f.data(), dy_f.data(),
                 wa_f.data(), Wa_f.data(),
                 R_invvar_f.data(), nanal_, neig_,
                 getkf_inflation, denkf, getkf, infl);
  this->Wa_ = Wa_f.cast<double>();
  this->wa_ = wa_f.cast<double>();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GETKFSolver<MODEL, OBS>::applyWeights(const IncrementEnsemble4D_ & bkg_pert,
                                           IncrementEnsemble4D_ & ana_pert,
                                           const GeometryIterator_ & i) {
  // apply Wa_, wa_
  util::Timer timer(classname(), "applyWeights");

  // allocate tmp arrays
  Eigen::MatrixXd XbModulated;  // modulated perturbations
  Eigen::MatrixXd XbOriginal;   // original perturbations

  // loop through analysis times and ens. members
  for (unsigned itime=0; itime < bkg_pert[0].size(); ++itime) {
    // cast bkg_pert ensemble at grid point i as an Eigen matrix Xb
    // modulates Xb
    XbModulated = vertloc_.modulateIncrement(bkg_pert, i, itime);
    // original Xb
    bkg_pert.packEigen(XbOriginal, i, itime);

    // postmulptiply
    // ensemble mean update
    Eigen::VectorXd xa = XbModulated*wa_;
    // ensemble perturbation update
    // Eq (10) from Lei 2018. (-) sign is accounted for in the Wa_ computation
    Eigen::MatrixXd Xa = XbOriginal + XbModulated*Wa_;

    // posterior inflation if rtps and rttp coefficients belong to (0,1]
    this->posteriorInflation(XbOriginal, Xa);

    // assign Xa_ to ana_pert
    Xa = Xa.colwise() + xa;
    ana_pert.setEigen(Xa, i, itime);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GETKFSolver<MODEL, OBS>::measurementUpdate(const IncrementEnsemble4D_ & bkg_pert,
                                                const GeometryIterator_ & i,
                                                IncrementEnsemble4D_ & ana_pert) {
  util::Timer timer(classname(), "measurementUpdate");

  // create the local subset of observations
  Departures_ locvector(this->obspaces_);
  locvector.ones();
  this->obsloc().computeLocalization(i, locvector);
  for (size_t iens = 0; iens < nanal_; ++iens) {
     (this->invVarR_)->mask(this->HZb_[iens]);
  }
  locvector.mask(*(this->invVarR_));
  Eigen::VectorXd local_omb_vec = this->omb_.packEigen(locvector);

  if (local_omb_vec.size() == 0) {
    // no obs. so no need to update Wa_ and wa_
    // ana_pert[i]=bkg_pert[i]
    this->copyLocalIncrement(bkg_pert, i, ana_pert);
  } else {
    // if obs are present do normal KF update
    // get local Yb & HZ
    Eigen::MatrixXd local_Yb_mat = this->Yb_.packEigen(locvector);
    Eigen::MatrixXd local_HZ_mat = this->HZb_.packEigen(locvector);
    // create local obs errors
    Eigen::VectorXd local_invVarR_vec = this->invVarR_->packEigen(locvector);
    // and apply localization
    Eigen::VectorXd localization = locvector.packEigen(locvector);
    local_invVarR_vec.array() *= localization.array();
    computeWeights(local_omb_vec, local_HZ_mat, local_Yb_mat, local_invVarR_vec);
    applyWeights(bkg_pert, ana_pert, i);
  }
}

}  // namespace oops
#endif  // OOPS_ASSIMILATION_GETKFSOLVER_H_
