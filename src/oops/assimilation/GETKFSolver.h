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
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/State4D.h"
#include "oops/base/StateEnsemble4D.h"
#include "oops/generic/VerticalLocEV.h"
#include "oops/interface/GeometryIterator.h"
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
  typedef IncrementEnsemble4D<MODEL>  IncrementEnsemble4D_;
  typedef ObsEnsemble<OBS>            ObsEnsemble_;
  typedef Observations<OBS>           Observations_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
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

  // compute/read H(x) for the original ensemble members
  // also computes omb_
  Observations_ yb_mean =
            LocalEnsembleSolver<MODEL, OBS>::computeHofX(ens_xx, iteration, readFromFile);
  if (readFromFile) {
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
    // modulate ensemble of obs
    IncrementEnsemble4D_ dx(ens_xx, this->xbmean_, this->incvars_);
    IncrementEnsemble4D_ Ztmp(geometry_, this->incvars_,
                              ens_xx[0].validTimes(), neig_);
    size_t ii = 0;
    for (size_t iens = 0; iens < nens_; ++iens) {
      Log::info() << " GETKFSolver::computeHofX starting ensemble member " << iens+1 << std::endl;
      util::printRunStats("GETKFSolver calculate hofx");
      vertloc_.modulateIncrement(dx[iens], Ztmp);
      for (size_t ieig = 0; ieig < neig_; ++ieig) {
        State4D_ tmpState = this->xbmean_;
        tmpState += Ztmp[ieig];
        Observations_ tmpObs(this->obspaces_);
        eckit::LocalConfiguration config;
        config.set("save hofx", false);
        config.set("save qc", false);
        config.set("save obs errors", false);
        this->computeHofX4D(config, tmpState, tmpObs);
        HZb_[ii] = tmpObs - yb_mean;
        tmpObs.save("hofxm"+std::to_string(iteration)+"_"+std::to_string(ieig+1)+
                      "_"+std::to_string(iens+1));
        ii = ii + 1;
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
