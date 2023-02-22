/*
 * (C) Copyright 2020 UCAR.
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

#include "oops/assimilation/LocalEnsembleSolver.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/Geometry.h"
#include "oops/base/IncrementEnsemble4D.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/ObsLocalizations.h"
#include "oops/base/ObsSpaces.h"
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
class LETKFSolver : public LocalEnsembleSolver<MODEL, OBS> {
  typedef Departures<OBS>             Departures_;
  typedef DeparturesEnsemble<OBS>     DeparturesEnsemble_;
  typedef Geometry<MODEL>             Geometry_;
  typedef GeometryIterator<MODEL>     GeometryIterator_;
  typedef IncrementEnsemble4D<MODEL>  IncrementEnsemble4D_;
  typedef ObsErrors<OBS>              ObsErrors_;
  typedef ObsLocalizations<MODEL, OBS> ObsLocalizations_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef State4D<MODEL>              State4D_;

 public:
  static const std::string classname() {return "oops::LETKFSolver";}

  LETKFSolver(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t,
              const State4D_ &, const Variables &);

  /// KF update + posterior inflation at a grid point location (GeometryIterator_)
  void measurementUpdate(const IncrementEnsemble4D_ &,
                         const GeometryIterator_ &, IncrementEnsemble4D_ &) override;

 protected:
  /// Computes weights for ensemble update with local observations
  /// \param[in] omb      Observation departures (nlocalobs)
  /// \param[in] Yb       Ensemble perturbations (nens, nlocalobs)
  /// \param[in] invvarR  Inverse of observation error variances (nlocalobs)
  virtual void computeWeights(const Eigen::VectorXd & omb, const Eigen::MatrixXd & Yb,
                              const Eigen::VectorXd & invvarR);

  /// Applies weights and adds posterior inflation
  virtual void applyWeights(const IncrementEnsemble4D_ &, IncrementEnsemble4D_ &,
                            const GeometryIterator_ &);

  Eigen::MatrixXd Wa_;  // transformation matrix for ens. perts. Xa=Xf*Wa
  Eigen::VectorXd wa_;  // transformation matrix for ens. mean xa=xf*wa

  // eigen solver matrices
  Eigen::VectorXd eival_;
  Eigen::MatrixXd eivec_;

  const size_t nens_;   // ensemble size
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
LETKFSolver<MODEL, OBS>::LETKFSolver(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                     const eckit::Configuration & config, size_t nens,
                                     const State4D_ & xbmean, const Variables & incvars)
  : LocalEnsembleSolver<MODEL, OBS>(obspaces, geometry, config, nens, xbmean, incvars),
    nens_(nens)
{
  Log::trace() << "LETKFSolver<MODEL, OBS>::create starting" << std::endl;
  Log::info() << "Using EIGEN implementation of LETKF" << std::endl;

  // pre-allocate transformation matrices
  Wa_.resize(nens_, nens_);
  wa_.resize(nens_);

  // pre-allocate eigen sovler matrices
  eival_.resize(nens_);
  eivec_.resize(nens_, nens_);
  Log::trace() << "LETKFSolver<MODEL, OBS>::create done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LETKFSolver<MODEL, OBS>::measurementUpdate(const IncrementEnsemble4D_ & bkg_pert,
                                                const GeometryIterator_ & i,
                                                IncrementEnsemble4D_ & ana_pert) {
  util::Timer timer(classname(), "measurementUpdate");

  // create the local subset of observations
  Departures_ locvector(this->obspaces_);
  locvector.ones();
  this->obsloc().computeLocalization(i, locvector);
  locvector.mask(*(this->invVarR_));
  Eigen::VectorXd local_omb_vec = this->omb_.packEigen(locvector);

  if (local_omb_vec.size() == 0) {
    // no obs. so no need to update Wa_ and wa_
    // ana_pert[i]=bkg_pert[i]
    this->copyLocalIncrement(bkg_pert, i, ana_pert);
  } else {
    // if obs are present do normal KF update
    // create local Yb
    Eigen::MatrixXd local_Yb_mat = this->Yb_.packEigen(locvector);
    // create local obs errors
    Eigen::VectorXd local_invVarR_vec = this->invVarR_->packEigen(locvector);
    // and apply localization
    Eigen::VectorXd localization = locvector.packEigen(locvector);
    local_invVarR_vec.array() *= localization.array();
    computeWeights(local_omb_vec, local_Yb_mat, local_invVarR_vec);
    applyWeights(bkg_pert, ana_pert, i);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LETKFSolver<MODEL, OBS>::computeWeights(const Eigen::VectorXd & dy,
                                             const Eigen::MatrixXd & Yb,
                                             const Eigen::VectorXd & diagInvR ) {
  // compute transformation matrix, save in Wa_, wa_
  // uses C++ eigen interface
  // implements LETKF from Hunt et al. 2007
  util::Timer timer(classname(), "computeWeights");

  const LocalEnsembleSolverInflationParameters & inflopt = this->options_.infl;

  // fill in the work matrix
  // work = Y^T R^-1 Y + (nens-1)/infl I
  const double infl = inflopt.mult;
  Eigen::MatrixXd work = Yb*(diagInvR.asDiagonal()*Yb.transpose());
  work.diagonal() += Eigen::VectorXd::Constant(nens_, (nens_-1)/infl);

  // eigenvalues and eigenvectors of the above matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(work);
  eival_ = es.eigenvalues().real();
  eivec_ = es.eigenvectors().real();

  // Pa   = [ Yb^T R^-1 Yb + (nens-1)/infl I ] ^-1
  work = eivec_ * eival_.cwiseInverse().asDiagonal() * eivec_.transpose();

  // Wa = sqrt[ (nens-1) Pa ]
  Wa_ = eivec_
      * ((nens_-1) * eival_.array().inverse()).sqrt().matrix().asDiagonal()
      * eivec_.transpose();

  // wa = Pa Yb^T R^-1 dy
  wa_ = work * (Yb * (diagInvR.asDiagonal()*dy));
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LETKFSolver<MODEL, OBS>::applyWeights(const IncrementEnsemble4D_ & bkg_pert,
                                           IncrementEnsemble4D_ & ana_pert,
                                           const GeometryIterator_ & i) {
  // applies Wa_, wa_
  util::Timer timer(classname(), "applyWeights");

  // loop through analysis times and ens. members
  for (size_t itime=0; itime < bkg_pert[0].size(); ++itime) {
    // make grid point forecast pert ensemble array
    Eigen::MatrixXd Xb;
    bkg_pert.packEigen(Xb, i, itime);

    // postmulptiply
    Eigen::VectorXd xa = Xb*wa_;   // ensemble mean update
    Eigen::MatrixXd Xa = Xb*Wa_;   // ensemble perturbation update

    // posterior inflation if rtps and rttp coefficients belong to (0,1]
    this->posteriorInflation(Xb, Xa);

    // assign Xa to ana_pert
    Xa = Xa.colwise() + xa;
    ana_pert.setEigen(Xa, i, itime);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LETKFSOLVER_H_
