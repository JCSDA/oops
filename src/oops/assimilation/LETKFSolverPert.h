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
class LETKFSolverPert : public LETKFSolver<MODEL, OBS> {
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
  static const std::string classname() {return "oops::LETKFSolverPert";}

  /// Constructor (instantiate LETKFSolver)
  LETKFSolverPert(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t,
              const StateSet_ &, const Variables &);

  Observations_ computeHofX(const StateEnsemble4D_ &, size_t, bool) override;

  /// KF update + posterior inflation at a grid point location (GeometryIterator_)
  void measurementUpdate(const IncrementEnsemble4D_ &,
                         const GeometryIterator_ &, IncrementEnsemble4D_ &) override;

 private:
  /// Computes weights for ensemble update with local observations
  /// \param[in] omb      Observation minus ensemble hofx mean (nlocalobs)
  /// \param[in] OmbPert  Observation perturbations minus ensemble hofx perturbations
  ///                     (nens, nlocalobs)
  /// \param[in] Yb       Ensemble hofx perturbations (nens, nlocalobs)
  /// \param[in] invVarR  Inverse of observation error variances (nlocalobs)
  virtual void computeWeights(const Eigen::VectorXd & omb, const Eigen::MatrixXd & OmbPert,
                              const Eigen::MatrixXd & Yb, const Eigen::VectorXd & invVarR);

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
LETKFSolverPert<MODEL, OBS>::LETKFSolverPert(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                     const eckit::Configuration & config, size_t nens,
                                     const StateSet_ & xbmean, const Variables & incvars)
  : LETKFSolver<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars),
    OmbPertDepEns_(obspaces, this->nens_)
{
  Log::trace() << "LETKFSolverPert<MODEL, OBS>::create starting" << std::endl;
  Log::trace() << "LETKFSolverPert<MODEL, OBS>::create done" << std::endl;
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
Observations<OBS> LETKFSolverPert<MODEL, OBS>::computeHofX(const StateEnsemble4D_ & ens_xx,
                                                       size_t iteration, bool readFromFile) {
  util::Timer timer(classname(), "computeHofX");

  Observations_ yb_mean(this->obspaces_);
  yb_mean = LETKFSolver<MODEL, OBS>::computeHofX(ens_xx, iteration, readFromFile);

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
void LETKFSolverPert<MODEL, OBS>::measurementUpdate(const IncrementEnsemble4D_ & bkg_pert,
                                                const GeometryIterator_ & i,
                                                IncrementEnsemble4D_ & ana_pert) {
  util::Timer timer(classname(), "measurementUpdate");

  // create the local subset of observations
  Departures_ locvector(this->obspaces_);
  locvector.ones();
  this->obsloc().computeLocalization(i, locvector);
  this->applyAssimilatedMask(locvector);
  const Eigen::VectorXd local_omb_vec = this->omb_.packEigen(locvector);

  if (local_omb_vec.size() == 0) {
    // no obs. so no need to update Wa_ and wa_
    // ana_pert[i]=bkg_pert[i]
    this->copyLocalIncrement(bkg_pert, i, ana_pert);
  } else {
    // if obs are present do normal KF update
    // create local OmbPert, Yb
    Eigen::MatrixXf local_OmbPert_mat_f = OmbPertDepEns_.packEigen(locvector);
    Eigen::MatrixXf local_Yb_mat_f = this->Yb_.packEigen(locvector);
    const Eigen::MatrixXd local_OmbPert_mat = local_OmbPert_mat_f.cast<double>();
    const Eigen::MatrixXd local_Yb_mat = local_Yb_mat_f.cast<double>();
    // create local obs errors and apply localization
    const Eigen::VectorXd localization = locvector.packEigen(locvector);
    const Eigen::VectorXd local_invVarR_vec = this->invVarR_->packEigen(locvector).array()
                                              * localization.array();
    computeWeights(local_omb_vec, local_OmbPert_mat, local_Yb_mat, local_invVarR_vec);
    applyWeights(bkg_pert, ana_pert, i);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LETKFSolverPert<MODEL, OBS>::computeWeights(const Eigen::VectorXd & omb,
                                                 const Eigen::MatrixXd & OmbPert,
                                                 const Eigen::MatrixXd & Yb,
                                                 const Eigen::VectorXd & invVarR ) {
  // compute transformation matrix, save in Wa_
  // uses C++ eigen interface
  // implements perturbed observation version of LETKF from Hunt et al. 2007
  util::Timer timer(classname(), "computeWeights");

  const LocalEnsembleSolverInflationParameters & inflopt = this->options_.infl;
  const double infl = inflopt.mult;

  // fill in the work matrix
  // work = Y^T R^-1 Y + (nens-1)/infl I
  Eigen::MatrixXd work = Yb*(invVarR.asDiagonal()*Yb.transpose());
  work.diagonal() += Eigen::VectorXd::Constant(this->nens_, ((this->nens_)-1)/infl);

  // eigenvalues and eigenvectors of the above matrix
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(work);
  this->eival_ = es.eigenvalues().real();
  this->eivec_ = es.eigenvectors().real();

  // Pa   = [ Yb^T R^-1 Yb + (nens-1)/infl I ] ^-1
  work.noalias() = (this->eivec_) * (this->eival_).cwiseInverse().asDiagonal()
                   * (this->eivec_).transpose();

  // Wa = Pa Yb^T R^-1 (dyPert)
  // dyPert = (y_mean + y_pert) - (yb_mean + yb_pert)
  (this->Wa_).noalias() = work * (Yb * (invVarR.asDiagonal()
                          * (OmbPert.transpose().colwise() + omb)));
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LETKFSolverPert<MODEL, OBS>::applyWeights(const IncrementEnsemble4D_ & bkg_pert,
                                           IncrementEnsemble4D_ & ana_pert,
                                           const GeometryIterator_ & i) {
  // applies Wa_, wa_
  util::Timer timer(classname(), "applyWeights");

  // loop through analysis times and ens. members
  for (size_t itime=0; itime < bkg_pert[0].size(); ++itime) {
    // make grid point forecast pert ensemble array
    Eigen::MatrixXd Xb;
    bkg_pert.packEigen(Xb, i, itime);

    // increment for ensemble members
    const Eigen::MatrixXd Xinc = Xb*(this->Wa_);
    // increment for ensemble mean
    const Eigen::VectorXd xinc = Xinc.rowwise().mean();
    // generate analysis perturbations for inflation
    Eigen::MatrixXd Xa = (Xinc + Xb).colwise() - xinc;

    // posterior inflation if rtps and rttp coefficients belong to (0,1]
    this->posteriorInflation(Xb, Xa);

    // add xinc to Xa and assign Xa to ana_pert
    Xa.colwise() += xinc;
    ana_pert.setEigen(Xa, i, itime);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LETKFSOLVERPERT_H_
