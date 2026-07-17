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
#include <numeric>
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
#include "oops/base/IncrementSet.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsLocalizations.h"
#include "oops/base/ObsSpaces.h"
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
  typedef IncrementSet<MODEL>         IncrementSet_;
  typedef LinearModel<MODEL>          LinearModel_;
  typedef Model<MODEL>                Model_;
  typedef ModelAuxControl<MODEL>      ModelAux_;
  typedef ModelAuxIncrement<MODEL>    ModelAuxInc_;
  typedef ObsAuxControls<OBS>         ObsAux_;
  typedef ObsAuxIncrements<OBS>       ObsAuxInc_;
  typedef ObsDataVector<OBS, int>     ObsDataInt_;
  typedef ObsEnsemble<OBS>            ObsEnsemble_;
  typedef ObsError<OBS>               ObsError_;
  typedef ObsErrors<OBS>              ObsErrors_;
  typedef Observations<OBS>           Observations_;
  typedef Observers<MODEL, OBS>       Observers_;
  typedef ObserversTLAD<MODEL, OBS>   ObserversTLAD_;
  typedef ObsLocalizations<MODEL, OBS> ObsLocalizations_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef PseudoModelState4D<MODEL>   PseudoModel_;
  typedef PseudoLinearModelIncrement4D<MODEL> PseudoLinearModel_;
  typedef State<MODEL>                State_;
  typedef StateSet<MODEL>             StateSet_;
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

  /// entire KF update (computeWeights+applyWeights) for a grid point GeometryIterator_
  void measurementUpdate(const Eigen::VectorXd &,
                         const ObsErrors_ &,
                         const Departures_ &,
                         const IncrementSet_ &,
                         const GeometryIterator_ &,
                         IncrementSet_ &) override;

 protected:
  Eigen::MatrixXf Wa_;  // transformation matrix for ens. perts. Xa_=Xf*Wa
  Eigen::VectorXf wa_;  // transformation matrix for ens. mean xa_=xf*wa
  size_t nens_;
  const Geometry_ & geometry_;
  VerticalLocEV_ vertloc_;
  size_t neig_;
  size_t nanal_;
  bool fortranETKF_;
  const bool useSVD_;

  std::unique_ptr<DeparturesEnsemble_> HZb_;

 private:
  /// Computes weights for ensemble update with local observations
  /// \param[in] omb      Observation departures (nlocalobs)
  /// \param[in] Yb       Ensemble perturbations for all the background memebers
  ///                     (nens*neig, nlocalobs)
  /// \param[in] YbOrig   Ensemble perturbations for the members to be updated (nens, nlocalobs)
  /// \param[in] invVarR  Inverse of observation error variances (nlocalobs)
  virtual void computeWeights(const Eigen::VectorXd & omb,
                              const Eigen::MatrixXf & Yb,
                              const Eigen::MatrixXf & YbOrig,
                              const ObsErrors_ & R);

  /// Applies weights and adds posterior inflation
  virtual void applyWeights(const IncrementSet_ &,
                            IncrementSet_ &,
                            const GeometryIterator_ &);
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
DeterministicGETKF<MODEL, OBS>::DeterministicGETKF(ObsSpaces_ & obspaces,
                                                   const Geometry_ & geometry,
                                                   const eckit::Configuration & config, size_t nens,
                                                   const StateSet_ & xbmean,
                                                   const Variables & incvars)
  : LocalEnsembleSolver<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars),
    nens_(nens),
    geometry_(geometry),
    vertloc_(config.getSubConfiguration("local ensemble DA.vertical localization"), xbmean[0],
             incvars),
    neig_(vertloc_.neig()),
    nanal_(neig_ * nens_),
    fortranETKF_(config.getBool("local ensemble DA.fortran ETKF", true)),
    useSVD_(this->svdRequested(config)) {
  Log::trace() << "DeterministicGETKF<MODEL, OBS>::create starting" << std::endl;
  // pre-allocate transformation matrices
  Wa_.resize(nanal_, nens);
  wa_.resize(nanal_);

  const eckit::LocalConfiguration localEnsConfig = config.getSubConfiguration("local ensemble DA");
  this->doCrossValidation = localEnsConfig.has("cross validation");

  if (this->doCrossValidation) {
    Log::trace() << "DeterministicGETKF<MODEL, OBS>::constructSubensembleSplitter starting"
                 << std::endl;
    const eckit::LocalConfiguration crossValidationConfig = localEnsConfig.getSubConfiguration(
                                                                           "cross validation");
    this->SubensembleSplitter_ = std::make_unique<oops::SubensembleSplitter>(this->nens_,
                                 this->neig_, crossValidationConfig);
    this->SubensembleSplitter_->split();
    this->nsubens_ = crossValidationConfig.getUnsigned("number of subensembles");
    Log::trace() << "DeterministicGETKF<MODEL, OBS>::constructSubensembleSplitter done"
                 << std::endl;
  }
  // initialize HZb_
  HZb_ = std::make_unique<DeparturesEnsemble_>(this->obspaces_, nanal_);

  // read modulated ensemble
  Observations_ ytmp(this->ybmean_);
  size_t ii = 0;
  for (size_t iens = 0; iens < nens_; ++iens) {
    Log::info() << " DeterministicGETKF::computeHofX starting ensemble member "
                << iens+1 << std::endl;
    util::printRunStats("DeterministicGETKF read hofx");
    for (size_t ieig = 0; ieig < neig_; ++ieig) {
      ytmp.read("hofxm0_"+std::to_string(ieig+1)+
                    "_"+std::to_string(iens+1));
      HZb_->setData(ii, ytmp - this->ybmean_);
      ii = ii + 1;
    }
  }

  Log::trace() << "DeterministicGETKF<MODEL, OBS>::create done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void DeterministicGETKF<MODEL, OBS>::computeWeights(const Eigen::VectorXd & dy,
                                                    const Eigen::MatrixXf & Yb,
                                                    const Eigen::MatrixXf & YbOrig,
                                                    const ObsErrors_ & R) {
  // compute transformation matrix, save in Wa_, wa_
  // Yb(nobs,neig*nens), YbOrig(nobs,nens)
  util::Timer timer(classname(), "computeWeights");
  const float infl = this->inflopt_.getFloat("mult", 1.0);

  if (fortranETKF_) {
    // cast eigen<double> to eigen<float>
    const Eigen::VectorXf dy_f = dy.cast<float>();
    const Eigen::VectorXd invVarR = R.local_invVarR();
    const Eigen::VectorXf invVarR_f = invVarR.cast<float>();

    // call into GSI interface to compute Wa and wa
    const int nobsl = dy.size();
    const int getkf_inflation = 0;
    const int denkf = 0;
    const int getkf = 1;
    letkf_core_f90(nobsl, Yb.data(), YbOrig.data(), dy_f.data(),
                   this->wa_.data(), this->Wa_.data(),
                   invVarR_f.data(), nanal_, neig_,
                   getkf_inflation, denkf, getkf, infl);
  } else {
    oops::detGETKF_computeWeights(dy.cast<float>(), Yb, YbOrig,
                                  R, infl, useSVD_, wa_, Wa_);
  }
}


// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void DeterministicGETKF<MODEL, OBS>::applyWeights(const IncrementSet_ & bkg_pert,
                                                  IncrementSet_ & ana_pert,
                                                  const GeometryIterator_ & i) {
  util::Timer timer(classname(), "applyWeights");

  oops::detETKF_applyWeights<MODEL>(bkg_pert,
                                    ana_pert,
                                    i,
                                    this->wa_,
                                    this->Wa_,
                                    this->inflopt_,
                                    &(vertloc_));
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void DeterministicGETKF<MODEL, OBS>::measurementUpdate(const Eigen::VectorXd & local_omb_vec,
                                                       const ObsErrors_ & R,
                                                       const Departures_ & locvector,
                                                       const IncrementSet_ & bkg_pert,
                                                       const GeometryIterator_ & i,
                                                       IncrementSet_ & ana_pert) {
  /*
  if (this->doCrossValidation) {
    throw eckit::NotImplemented(
          "Cross validation for deterministic GETKF not yet implemented",
          Here());
  }
  */
  const Eigen::MatrixXf local_Yb_mat_f = (this->Yb_)->packEigen(locvector);
  const Eigen::MatrixXf local_HZ_mat_f = (this->HZb_)->packEigen(locvector);
  this->computeWeights(local_omb_vec, local_HZ_mat_f, local_Yb_mat_f, R);
  this->applyWeights(bkg_pert, ana_pert, i);
}


}  // namespace oops
#endif  // OOPS_ASSIMILATION_GETKFSOLVER_H_
