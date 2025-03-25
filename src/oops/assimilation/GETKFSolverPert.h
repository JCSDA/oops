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
class GETKFSolverPert : public GETKFSolver<MODEL, OBS> {
  typedef Departures<OBS>             Departures_;
  typedef DeparturesEnsemble<OBS>     DeparturesEnsemble_;
  typedef Geometry<MODEL>             Geometry_;
  typedef GeometryIterator<MODEL>     GeometryIterator_;
  typedef Observations<OBS>           Observations_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef StateSet<MODEL>             StateSet_;
  typedef StateEnsemble4D<MODEL>      StateEnsemble4D_;
  typedef IncrementEnsemble4D<MODEL>  IncrementEnsemble4D_;

 public:
  static const std::string classname() {return "oops::GETKFSolverPert";}
  /// Constructor (instantiate GETKFSolver, eival and eivec)
  GETKFSolverPert(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t,
              const StateSet_ &, const Variables &);

  Observations_ computeHofX(const StateEnsemble4D_ &, size_t, bool) override;

  /// entire KF update (computeWeights+applyWeights) for a grid point GeometryIterator_
  void measurementUpdate(const IncrementEnsemble4D_ &, const GeometryIterator_ &,
                         IncrementEnsemble4D_ &) override;

 private:
  /// Computes weights for ensemble update with local observations
  /// \param[in] omb   Observation minus ensemble hofx mean (nlocalobs)
  /// \param[in] Yb    Observation perturbations minus original ensemble hofx perturbations
  ///                  (nens, nlocalobs)
  /// \param[in] HZb   Modulated ensemble hofx perturbations (nens*neig, nlocalobs)
  /// \param[in] invVarR  Inverse of observation error variances (nlocalobs)
  void computeWeights(const Eigen::VectorXd & omb, const Eigen::MatrixXd & OmbPert,
                      const Eigen::MatrixXd & HZb, const Eigen::VectorXd & invVarR);

  /// Applies weights and adds posterior inflation
  void applyWeights(const IncrementEnsemble4D_ &, IncrementEnsemble4D_ &,
                    const GeometryIterator_ &);

 private:
  // parameters
  Eigen::VectorXf eival_;
  Eigen::MatrixXf eivec_;
};

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
GETKFSolverPert<MODEL, OBS>::GETKFSolverPert(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                const eckit::Configuration & config, size_t nens,
                                const StateSet_ & xbmean, const Variables & incvars)
  : GETKFSolver<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars),
    eival_(this->nanal_), eivec_(this->nanal_, this->nanal_)
{
  Log::trace() << "GETKFSolverPert<MODEL, OBS>::create starting" << std::endl;
  Log::trace() << "GETKFSolverPert<MODEL, OBS>::create done" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
Observations<OBS> GETKFSolverPert<MODEL, OBS>::computeHofX(const StateEnsemble4D_ & ens_xx,
                                                       size_t iteration, bool readFromFile) {
  util::Timer timer(classname(), "computeHofX");

  Observations_ yb_mean(this->obspaces_);
  yb_mean = GETKFSolver<MODEL, OBS>::computeHofX(ens_xx, iteration, readFromFile);

  // generate observation perturbations with zero mean and
  // store in Yb_ observation perturbations minus original ensemble hofx
  // perturbations to be consistent with omb used in Kalman gain calculation
  Departures_ ypertDepTmp(this->obspaces_);
  Departures_ ypertDepSum(this->obspaces_);

  ypertDepSum.zero();
  for (size_t iens = 0; iens < (this->nens_); ++iens) {
    ypertDepTmp.zero();
    (*(this->R_)).randomize(ypertDepTmp);
    (this->Yb_)[iens] *= -1;
    (this->Yb_)[iens] += ypertDepTmp;
    ypertDepSum += ypertDepTmp;
  }

  for (size_t iens = 0; iens < (this->nens_); ++iens) {
    (this->Yb_)[iens].axpy(-1.0/(this->nens_), ypertDepSum);
  }

  return yb_mean;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GETKFSolverPert<MODEL, OBS>::computeWeights(const Eigen::VectorXd & omb,
                                                 const Eigen::MatrixXd & OmbPert,
                                                 const Eigen::MatrixXd & HZb,
                                                 const Eigen::VectorXd & invVarR) {
  // compute transformation matrix, save in Wa_
  // uses C++ eigen interface
  util::Timer timer(classname(), "computeWeights");
  const LocalEnsembleSolverInflationParameters & inflopt = this->options_.infl;
  const float infl = inflopt.mult;

  const Eigen::VectorXf omb_f = omb.cast<float>();
  const Eigen::MatrixXf OmbPert_f = OmbPert.cast<float>();
  const Eigen::MatrixXf HZb_f = HZb.cast<float>();
  const Eigen::VectorXf invVarR_f = invVarR.cast<float>();

  // fill in the work matrix
  // work = Y^T R^-1 Y + (nens-1)/infl I
  Eigen::MatrixXf work = HZb_f*(invVarR_f.asDiagonal()*HZb_f.transpose());
  work.diagonal() += Eigen::VectorXf::Constant(this->nanal_, ((this->nens_)-1)/infl);

  // compute eigenvalues and eigenvectors of work matrix
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(work);
  eival_ = es.eigenvalues().real();
  eivec_ = es.eigenvectors().real();

  // Pa   = [ Yb^T R^-1 Yb + (nens-1)/infl I ] ^-1
  work.noalias() = eivec_ * eival_.cwiseInverse().asDiagonal() * eivec_.transpose();

  // Wa = Pa Yb^T R^-1 (dyPert)
  // dyPert = (y_mean + y_pert) - (yb_mean + yb_pert)
  Eigen::MatrixXf Wa_f(this->nanal_, this->nens_);
  Wa_f.noalias() = work * (HZb_f * (invVarR_f.asDiagonal()
                   * (OmbPert_f.transpose().colwise() + omb_f)));

  this->Wa_ = Wa_f.cast<double>();
}

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GETKFSolverPert<MODEL, OBS>::applyWeights(const IncrementEnsemble4D_ & bkg_pert,
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
    XbModulated = this->vertloc_.modulateIncrement(bkg_pert, i, itime);
    // original Xb
    bkg_pert.packEigen(XbOriginal, i, itime);

    // increment for ensemble members
    const Eigen::MatrixXd Xinc = XbModulated*(this->Wa_);
    // increment for ensemble mean
    const Eigen::VectorXd xinc = Xinc.rowwise().mean();
    // generate analysis perturbations for inflation
    Eigen::MatrixXd Xa = (Xinc + XbOriginal).colwise() - xinc;

    // posterior inflation if rtps and rttp coefficients belong to (0,1]
    this->posteriorInflation(XbOriginal, Xa);

    // add xinc to Xa and assign Xa to ana_pert
    Xa.colwise() += xinc;
    ana_pert.setEigen(Xa, i, itime);
  }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GETKFSolverPert<MODEL, OBS>::measurementUpdate(const IncrementEnsemble4D_ & bkg_pert,
                                                const GeometryIterator_ & i,
                                                IncrementEnsemble4D_ & ana_pert) {
  util::Timer timer(classname(), "measurementUpdate");

  // create the local subset of observations
  Departures_ locvector(this->obspaces_);
  locvector.ones();
  this->obsloc().computeLocalization(i, locvector);
  locvector.mask(*(this->invVarR_));
  const Eigen::VectorXd local_omb_vec = this->omb_.packEigen(locvector);

  if (local_omb_vec.size() == 0) {
    // no obs. so no need to update Wa_ and wa_
    // ana_pert[i]=bkg_pert[i]
    this->copyLocalIncrement(bkg_pert, i, ana_pert);
  } else {
    // if obs are present do normal KF update
    // get local Yb, HZb
    // (this->Yb_) stores obs. perturbations minus
    // original ensemble hofx perturbations
    const Eigen::MatrixXd local_OmbPert_mat = (this->Yb_).packEigen(locvector);
    const Eigen::MatrixXd local_HZb_mat = (this->HZb_).packEigen(locvector);
    // create local obs errors and apply localization
    const Eigen::VectorXd localization = locvector.packEigen(locvector);
    const Eigen::VectorXd local_invVarR_vec = (this->invVarR_)->packEigen(locvector).array()
                                              * localization.array();
    computeWeights(local_omb_vec, local_OmbPert_mat, local_HZb_mat, local_invVarR_vec);
    applyWeights(bkg_pert, ana_pert, i);
  }
}

}  // namespace oops
#endif  // OOPS_ASSIMILATION_GETKFSOLVERPERT_H_
