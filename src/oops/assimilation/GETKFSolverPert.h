/*
 * (C) Copyright 2024-2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_GETKFSOLVERPERT_H_
#define OOPS_ASSIMILATION_GETKFSOLVERPERT_H_

#include <Eigen/Dense>
#include <cfloat>
#include <memory>
#include <string>
#include <vector>

namespace oops {
  class Variables;


/* An implementation of stochastic GETKF of assimilating
 * perturbed observations as a derived class of deterministic
 * GETKF solver (GETKFSolver.h)
 */

template <typename MODEL, typename OBS>
class StochasticGETKF : public DeterministicGETKF<MODEL, OBS> {
  typedef Departures<OBS>             Departures_;
  typedef DeparturesEnsemble<OBS>     DeparturesEnsemble_;
  typedef Geometry<MODEL>             Geometry_;
  typedef GeometryIterator<MODEL>     GeometryIterator_;
  typedef Observations<OBS>           Observations_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef StateSet<MODEL>             StateSet_;
  typedef StateEnsemble4D<MODEL>      StateEnsemble4D_;
  typedef IncrementSet<MODEL>         IncrementSet_;
  typedef VerticalLocEV<MODEL>        VerticalLocEV_;

 public:
  static const std::string classname() {return "oops::StochasticGETKF";}
  /// Constructor (instantiate GETKFSolver, eival and eivec)
  StochasticGETKF(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t,
                  const StateSet_ &, const Variables &);

  Observations_ computeHofX(const StateEnsemble4D_ &, size_t, bool) override;

  /// entire KF update (computeWeights+applyWeights) for a grid point GeometryIterator_
  void measurementUpdate(const Eigen::VectorXd &,
                         const Eigen::VectorXd &,
                         const Departures_ &,
                         const IncrementSet_ &,
                         const GeometryIterator_ &,
                         IncrementSet_ &) override;

 private:
  /// Computes weights for ensemble update with local observations
  /// \param[in] omb   Observation minus ensemble hofx mean (nlocalobs)
  /// \param[in] Yb    Observation perturbations minus original ensemble hofx perturbations
  ///                  (nens, nlocalobs)
  /// \param[in] HZb   Modulated ensemble hofx perturbations (nens*neig, nlocalobs)
  /// \param[in] invVarR  Inverse of observation error variances (nlocalobs)
  void computeWeights(const Eigen::VectorXd & omb,
                      const Eigen::MatrixXf & OmbPert_f,
                      const Eigen::MatrixXf & HZb_f,
                      const Eigen::VectorXd & invVarR) override;

  /// Applies weights and adds posterior inflation
  void applyWeights(const IncrementSet_ &,
                    IncrementSet_ &,
                    const GeometryIterator_ &) override;

 private:
  // parameters
  Eigen::VectorXf eival_;
  Eigen::MatrixXf eivec_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
StochasticGETKF<MODEL, OBS>::StochasticGETKF(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                             const eckit::Configuration & config, size_t nens,
                                             const StateSet_ & xbmean, const Variables & incvars)
  : DeterministicGETKF<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars),
    eival_(this->nanal_), eivec_(this->nanal_, this->nanal_)
{
  Log::trace() << "StochasticGETKF<MODEL, OBS>::create starting" << std::endl;
  Log::trace() << "StochasticGETKF<MODEL, OBS>::create done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observations<OBS> StochasticGETKF<MODEL, OBS>::computeHofX(const StateEnsemble4D_ & ens_xx,
                                                           size_t iteration, bool readFromFile) {
  util::Timer timer(classname(), "computeHofX");

  Observations_ yb_mean(this->obspaces_);
  yb_mean = DeterministicGETKF<MODEL, OBS>::computeHofX(ens_xx, iteration, readFromFile);

  // generate observation perturbations with zero mean and
  // store in Yb_ observation perturbations minus original ensemble hofx
  // perturbations to be consistent with omb used in Kalman gain calculation
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
    // Temporary storage in Yb_[i]
    (this->Yb_)->setData(iens, pertDepTmp);
  }

  for (size_t iens = 0; iens < (this->nens_); ++iens) {
    pertDepTmp.zero();
    pertDepTmp = (this->Yb_)->getData(iens);
    // Subtract mean(dy_j) from dy_i - Y_i
    pertDepTmp.axpy(-1.0/(this->nens_), ypertDepSum);
    // Now Yb_[i] = dy_i - Y_i - mean(dy_j) forall j in ens
    (this->Yb_)->setData(iens, pertDepTmp);
  }

  return yb_mean;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void StochasticGETKF<MODEL, OBS>::computeWeights(const Eigen::VectorXd & dy,
                                                 const Eigen::MatrixXf & YbOrig,
                                                 const Eigen::MatrixXf & Yb,
                                                 const Eigen::VectorXd & invVarR) {
  // compute transformation matrix, save in Wa_
  util::Timer timer(classname(), "computeWeights");
  const LocalEnsembleSolverInflationParameters & inflopt = this->options_.infl;
  const float infl = inflopt.mult;

  oops::stoETKF_computeWeights(dy, Yb, YbOrig, invVarR, infl, this->Wa_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void StochasticGETKF<MODEL, OBS>::applyWeights(const IncrementSet_ & bkg_pert,
                                               IncrementSet_ & ana_pert,
                                               const GeometryIterator_ & i) {
  util::Timer timer(classname(), "applyWeights");

  oops::stoETKF_applyWeights<MODEL>(bkg_pert,
                                    ana_pert,
                                    i,
                                    this->Wa_,
                                    this->options_.infl,
                                    &(this->vertloc_));
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void StochasticGETKF<MODEL, OBS>::measurementUpdate(const Eigen::VectorXd & local_omb_vec,
                                                    const Eigen::VectorXd & local_invVarR_vec,
                                                    const Departures_ & locvector,
                                                    const IncrementSet_ & bkg_pert,
                                                    const GeometryIterator_ & i,
                                                    IncrementSet_ & ana_pert) {
  const Eigen::MatrixXf local_OmbPert_mat_f = (this->Yb_)->packEigen(locvector);
  const Eigen::MatrixXf local_HZb_mat_f = (this->HZb_)->packEigen(locvector);
  this->computeWeights(local_omb_vec, local_OmbPert_mat_f, local_HZb_mat_f, local_invVarR_vec);
  this->applyWeights(bkg_pert, ana_pert, i);
}

}  // namespace oops
#endif  // OOPS_ASSIMILATION_GETKFSOLVERPERT_H_
