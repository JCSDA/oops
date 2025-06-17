/*
 * (C) Copyright 2024-2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#ifndef OOPS_ASSIMILATION_LETKFSOLVERPERT_H_
#define OOPS_ASSIMILATION_LETKFSOLVERPERT_H_

#include <Eigen/Dense>
#include <cfloat>
#include <memory>
#include <string>
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
  typedef IncrementEnsemble4D<MODEL>  IncrementEnsemble4D_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef StateSet<MODEL>             StateSet_;
  typedef StateEnsemble4D<MODEL>      StateEnsemble4D_;

 public:
  static const std::string classname() {return "oops::StochasticLETKF";}

  /// Constructor (instantiate LETKFSolver)
  StochasticLETKF(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t,
                  const StateSet_ &, const Variables &);

  Observations_ computeHofX(const StateEnsemble4D_ &, size_t, bool) override;

  /// KF update + posterior inflation at a grid point location (GeometryIterator_)
  void measurementUpdate(const Eigen::VectorXd &,
                         const Eigen::VectorXd &,
                         const Departures_ &,
                         const IncrementEnsemble4D_ &,
                         const GeometryIterator_ &,
                         IncrementEnsemble4D_ &) override;

 private:
  /// Computes weights for ensemble update with local observations
  /// \param[in] omb      Observation minus ensemble hofx mean (nlocalobs)
  /// \param[in] OmbPert  Observation perturbations minus ensemble hofx perturbations
  ///                     (nens, nlocalobs)
  /// \param[in] Yb       Ensemble hofx perturbations (nens, nlocalobs)
  /// \param[in] invVarR  Inverse of observation error variances (nlocalobs)
  virtual void computeWeights(const Eigen::VectorXd & omb,
                              const Eigen::MatrixXf & OmbPert,
                              const Eigen::MatrixXf & Yb,
                              const Eigen::VectorXd & invVarR);

  /// Applies weights and adds posterior inflation
  virtual void applyWeights(const IncrementEnsemble4D_ &, IncrementEnsemble4D_ &,
                            const GeometryIterator_ &);

 private:
  // departure ensemble object of observation perturbations
  // minus ensemble hofx perturbations
  DeparturesEnsemble_ OmbPertDepEns_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
StochasticLETKF<MODEL, OBS>::StochasticLETKF(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                     const eckit::Configuration & config, size_t nens,
                                     const StateSet_ & xbmean, const Variables & incvars)
  : DeterministicLETKF<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars),
    OmbPertDepEns_(obspaces, this->nens_)
{
  Log::trace() << "StochasticLETKF<MODEL, OBS>::create starting" << std::endl;
  Log::trace() << "StochasticLETKF<MODEL, OBS>::create done" << std::endl;
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
Observations<OBS> StochasticLETKF<MODEL, OBS>::computeHofX(const StateEnsemble4D_ & ens_xx,
                                                       size_t iteration, bool readFromFile) {
  util::Timer timer(classname(), "computeHofX");

  Observations_ yb_mean(this->obspaces_);
  yb_mean = DeterministicLETKF<MODEL, OBS>::computeHofX(ens_xx, iteration, readFromFile);

  // generate observation perturbations with zero mean and
  // compute observation perturbations minus ensemble hofx perturbations
  Departures_ pertDepTmp(this->obspaces_);
  Departures_ ypertDepSum(this->obspaces_);

  ypertDepSum.zero();
  for (size_t iens = 0; iens < (this->nens_); ++iens) {
    pertDepTmp.zero();
    (*(this->R_)).randomize(pertDepTmp);
    ypertDepSum += pertDepTmp;
    pertDepTmp -= (this->Yb_).getData(iens);
    OmbPertDepEns_.setData(iens, pertDepTmp);
  }

  for (size_t iens = 0; iens < (this->nens_); ++iens) {
    pertDepTmp.zero();
    pertDepTmp = OmbPertDepEns_.getData(iens);
    pertDepTmp.axpy(-1.0/(this->nens_), ypertDepSum);
    OmbPertDepEns_.setData(iens, pertDepTmp);
  }

  return yb_mean;
}


// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void StochasticLETKF<MODEL, OBS>::measurementUpdate(const Eigen::VectorXd & local_omb_vec,
                                                    const Eigen::VectorXd & local_invVarR_vec,
                                                    const Departures_ & locvector,
                                                    const IncrementEnsemble4D_ & bkg_pert,
                                                    const GeometryIterator_ & i,
                                                    IncrementEnsemble4D_ & ana_pert) {
  const Eigen::MatrixXf local_OmbPert_mat_f = OmbPertDepEns_.packEigen(locvector);
  const Eigen::MatrixXf local_Yb_mat_f = this->Yb_.packEigen(locvector);
  this->computeWeights(local_omb_vec, local_OmbPert_mat_f, local_Yb_mat_f, local_invVarR_vec);
  this->applyWeights(bkg_pert, ana_pert, i);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void StochasticLETKF<MODEL, OBS>::computeWeights(const Eigen::VectorXd & dy,
                                                 const Eigen::MatrixXf & YbOrig,
                                                 const Eigen::MatrixXf & Yb,
                                                 const Eigen::VectorXd & invVarR) {
  // compute transformation matrix, save in Wa_
  // implements perturbed observation version of LETKF from Hunt et al. 2007
  util::Timer timer(classname(), "computeWeights");
  const LocalEnsembleSolverInflationParameters & inflopt = this->options_.infl;
  const double infl = inflopt.mult;

  oops::stoETKF_computeWeights(dy, Yb, YbOrig, invVarR, infl, this->Wa_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void StochasticLETKF<MODEL, OBS>::applyWeights(const IncrementEnsemble4D_ & bkg_pert,
                                               IncrementEnsemble4D_ & ana_pert,
                                               const GeometryIterator_ & i) {
  util::Timer timer(classname(), "applyWeights");

  oops::stoETKF_applyWeights<MODEL>(bkg_pert,
                                    ana_pert,
                                    i,
                                    this->Wa_,
                                    this->options_.infl);
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LETKFSOLVERPERT_H_
