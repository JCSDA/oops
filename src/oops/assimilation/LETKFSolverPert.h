/*
 * (C) Copyright 2024-2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#ifndef OOPS_ASSIMILATION_LETKFSOLVERPERT_H_
#define OOPS_ASSIMILATION_LETKFSOLVERPERT_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cfloat>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace oops {
  class Variables;

/* An implementation of stochastic LETKF of assimilating
 * pertrubed observations as a derived class of deterministic
 * LETKF solver (LETKFSolver.h)
 */
template <typename MODEL, typename OBS>
class StochasticLETKF : public DeterministicLETKF<MODEL, OBS> {
  typedef Departures<OBS>             Departures_;
  typedef DeparturesEnsemble<OBS>     DeparturesEnsemble_;
  typedef Geometry<MODEL>             Geometry_;
  typedef GeometryIterator<MODEL>     GeometryIterator_;
  typedef Observations<OBS>           Observations_;
  typedef IncrementSet<MODEL>         IncrementSet_;
  typedef ObsErrors<OBS>              ObsErrors_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef StateSet<MODEL>             StateSet_;

 public:
  static const std::string classname() {return "oops::StochasticLETKF";}

  /// Constructor (instantiate LETKFSolver)
  StochasticLETKF(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t,
                  const StateSet_ &, const Variables &);

  /// KF update + posterior inflation at a grid point location (GeometryIterator_)
  void measurementUpdate(const Eigen::VectorXd &,
                         const ObsErrors_ &,
                         const Departures_ &,
                         const IncrementSet_ &,
                         const GeometryIterator_ &,
                         IncrementSet_ &) override;

 private:
  /// Computes weights for ensemble update with local observations
  /// \param[in] omb      Observation minus ensemble hofx mean (nlocalobs)
  /// \param[in] OmbPert  Observation perturbations minus ensemble hofx perturbations
  ///                     (nens, nlocalobs)
  /// \param[in] Yb       Ensemble hofx perturbations (nens, nlocalobs)
  /// \param[in] invVarR  Inverse of observation error variances (nlocalobs)
  void computeWeights(const Eigen::VectorXd & omb,
                      const Eigen::MatrixXf & OmbPert,
                      const Eigen::MatrixXf & Yb,
                      const ObsErrors_ & R);

  /// Computes weights for ensemble update with local observations for a given projection matrix
  /// \param[in] omb                 Observation minus ensemble hofx mean (nlocalobs)
  /// \param[in] OmbPert             Observation perturbations minus ensemble hofx perturbations
  ///                                (nens, nlocalobs)
  /// \param[in] YbRinvYbpI          Matrix equal to Y^T R^-1 Y + (nens-1)/infl I (nens, nens)
  /// \param[in] YbRinv              Observation perturbations minus original ensemble hofx
  ///                                perturbations multiplying inverse R (nens, nlocalobs)
  /// \param[in] excludedProjection  Projection matrix for the excluded
  ///                                subensemble members for cross validation (nens, nhat)
  /// \param[in] includedProjection  Projection matrix for the included
  ///                                subensemble members for cross validation (nens, nens)
  virtual void computeWeights(const Eigen::VectorXd & omb,
                              const Eigen::MatrixXf & OmbPert,
                              const Eigen::MatrixXf & YbRinvYbpI,
                              const Eigen::MatrixXf & YbRinv,
                              const Eigen::SparseMatrix<float> & excludedProjection,
                              const Eigen::SparseMatrix<float> & includedProjection);

  /// Computes localised YbRinv and YbRinvYbpI
  /// \param[in]  locvector          Departures vector used for localisation
  /// \param[in]  local_invVarR_vec  Localised inverse variance of the R matrix
  ///                                (nlocalobs, nlocalobs)
  /// \param[out] YbRinv             Observation perturbations minus original ensemble hofx
  ///                                perturbations multiplying inverse R (nens, nlocalobs)
  /// \param[out] YbRinvYbpI         Matrix equal to Y^T R^-1 Y + (nens-1)/infl I
  ///                                (nens, nens)
  const std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> computeYbRinvMatrices(const Departures_ &
                                                                           locvector,
                                                                           const ObsErrors_ & R);

  /// Applies weights and adds posterior inflation
  void applyWeights(const IncrementSet_ &, IncrementSet_ &,
                    const GeometryIterator_ &) override;

 private:
  // departure ensemble object of observation perturbations
  // minus ensemble hofx perturbations
  DeparturesEnsemble_ OmbPertDepEns_;
  const bool useSVD_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
StochasticLETKF<MODEL, OBS>::StochasticLETKF(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                             const eckit::Configuration & config, size_t nens,
                                             const StateSet_ & xbmean, const Variables & incvars)
  : DeterministicLETKF<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars),
    OmbPertDepEns_(obspaces, this->nens_),
    useSVD_(this->svdRequested(config)) {
  Log::trace() << "StochasticLETKF<MODEL, OBS>::create starting" << std::endl;
  // generate observation perturbations with zero mean and
  // compute observation perturbations minus ensemble hofx perturbations
  Departures_ pertDepTmp(this->obspaces_);
  Departures_ ypertDepSum(this->obspaces_);

  ypertDepSum.zero();
  for (size_t iens = 0; iens < (this->nens_); ++iens) {
    pertDepTmp.zero();
    if (static_cast<int>(iens) != this->unperturbedIdx_) {
      // Get perturbations dy_i (zero for unperturbed member)
      (*(this->R_)).randomize(pertDepTmp);
    }
    // Sum the perturbations for calculating the mean perturbation mean(dy_j)
    // forall j in ens later.
    ypertDepSum += pertDepTmp;
    // Now pertDepTmp = dy_i - Y_i
    //                = dy_i - (H(x_i) - mean(H(x_j))) forall j in ens
    pertDepTmp -= (this->Yb_)->getData(iens);
    // Temporary storage in OmbPertDepEns_[i]
    OmbPertDepEns_.setData(iens, pertDepTmp);
  }

  for (size_t iens = 0; iens < (this->nens_); ++iens) {
    pertDepTmp.zero();
    pertDepTmp = OmbPertDepEns_.getData(iens);
    // Subtract mean(dy_j) from dy_i - Y_i
    pertDepTmp.axpy(-1.0/(this->nens_), ypertDepSum);
    // Now OmbPertDepEns_[i] = dy_i - Y_i - mean(dy_j) forall j in ens
    OmbPertDepEns_.setData(iens, pertDepTmp);
  }
  Log::trace() << "StochasticLETKF<MODEL, OBS>::create done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void StochasticLETKF<MODEL, OBS>::measurementUpdate(const Eigen::VectorXd & local_omb_vec,
                                                const ObsErrors_ & R,
                                                const Departures_ & locvector,
                                                const IncrementSet_ & bkg_pert,
                                                const GeometryIterator_ & i,
                                                IncrementSet_ & ana_pert) {
    const Eigen::MatrixXf local_OmbPert_mat_f = OmbPertDepEns_.packEigen(locvector);
    if (this->doCrossValidation) {
        this->Wa_.setZero();
        const std::tuple<Eigen::MatrixXf, Eigen::MatrixXf>
            ETKFCoreMatrices = this->computeYbRinvMatrices(locvector, R);
        const Eigen::MatrixXf & YbRinv = std::get<0>(ETKFCoreMatrices);
        const Eigen::MatrixXf & YbRinvYbpI = std::get<1>(ETKFCoreMatrices);
        const bool modulated = false;

    for (size_t isubens = 0; isubens < (this->nsubens_); ++isubens) {
      const std::tuple<Eigen::SparseMatrix<float>, Eigen::SparseMatrix<float>, std::vector<size_t>>
      projectionMatrices = this->SubensembleSplitter_->getProjectionMatrices(isubens, modulated);
      const Eigen::SparseMatrix<float> & excludedProjection = std::get<0>(projectionMatrices);
      const Eigen::SparseMatrix<float> & includedProjection = std::get<1>(projectionMatrices);
      this->computeWeights(local_omb_vec, local_OmbPert_mat_f, YbRinvYbpI, YbRinv,
                           excludedProjection, includedProjection);
    }
  } else {
    const Eigen::MatrixXf local_Yb_mat_f = (this->Yb_)->packEigen(locvector);
    this->computeWeights(local_omb_vec, local_OmbPert_mat_f, local_Yb_mat_f, R);
  }
  this->applyWeights(bkg_pert, ana_pert, i);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void StochasticLETKF<MODEL, OBS>::computeWeights(const Eigen::VectorXd & dy,
                                                 const Eigen::MatrixXf & YbOrig,
                                                 const Eigen::MatrixXf & Yb,
                                                 const ObsErrors_ & R) {
    // compute transformation matrix, save in Wa_
    // implements perturbed observation version of LETKF from Hunt et al. 2007
    util::Timer timer(classname(), "computeWeights");
    const float infl = this->inflopt_.getFloat("mult", 1.0);

  oops::stoETKF_computeWeights(dy.cast<float>(), Yb, YbOrig,
                               R, infl, useSVD_, this->Wa_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void StochasticLETKF<MODEL, OBS>::computeWeights(const Eigen::VectorXd & dy,
                                                 const Eigen::MatrixXf & YbOrig,
                                                 const Eigen::MatrixXf & YbRinvYbpI,
                                                 const Eigen::MatrixXf & YbRinv,
                                                 const Eigen::SparseMatrix<float> &
                                                 excludedProjection,
                                                 const Eigen::SparseMatrix<float> &
                                                 includedProjection) {
  // compute transformation matrix, save in Wa_
  // implements perturbed observation version of LETKF from Hunt et al. 2007
  util::Timer timer(classname(), "computeWeights");

  oops::stoETKF_computeWeights(dy.cast<float>(), YbRinvYbpI, YbRinv, YbOrig,
                               excludedProjection, includedProjection, this->Wa_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
const std::tuple<Eigen::MatrixXf, Eigen::MatrixXf>
StochasticLETKF<MODEL, OBS>::computeYbRinvMatrices(const Departures_ & locvector,
                                                   const ObsErrors_ & R) {
  // Pre-calculate YbRinv and YbRinvYbpI
  const Eigen::MatrixXf local_Yb_mat_f = (this->Yb_)->packEigen(locvector);
  const float infl = this->inflopt_.getFloat("mult", 1.0);
  const float scale = (this->nens_ - 1) / infl;

  const Eigen::MatrixXf YbRinv = oops::ETKF_YbRinv(local_Yb_mat_f, R);
  const Eigen::MatrixXf YbRinvYbpI = oops::ETKF_YbRinvYbpI(local_Yb_mat_f, YbRinv, scale);

    return std::make_tuple(YbRinv, YbRinvYbpI);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void StochasticLETKF<MODEL, OBS>::applyWeights(const IncrementSet_ & bkg_pert,
                                               IncrementSet_ & ana_pert,
                                               const GeometryIterator_ & i) {
  util::Timer timer(classname(), "applyWeights");

  oops::stoETKF_applyWeights<MODEL>(bkg_pert,
                                    ana_pert,
                                    i,
                                    this->Wa_,
                                    this->inflopt_);
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LETKFSOLVERPERT_H_
