/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

/// Local Ensemble Tranform Kalman Filter solver
/*!
 * An initial implementation of the LETKF from Hunt et al. 2007
 * this version is implemented using pure OOPS primitives with no
 * no temporary Eigen matrices use for Xa and Xb
 * this verion implemnts RTPP but not RTPS
 *
 * Hunt, B. R., Kostelich, E. J., & Szunyogh, I. (2007). Efficient data
 * assimilation for spatiotemporal chaos: A local ensemble transform Kalman
 * filter. Physica D: Nonlinear Phenomena, 230(1-2), 112-126.
 */


#ifndef OOPS_ASSIMILATION_LETKFSOLVEROOPS_H_
#define OOPS_ASSIMILATION_LETKFSOLVEROOPS_H_

#include <Eigen/Dense>

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/LETKFSolver.h"
#include "oops/assimilation/LETKFSolverParameters.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/IncrementEnsemble4D.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {

template <typename MODEL, typename OBS>
class LETKFSolverOOPS : public LETKFSolver<MODEL, OBS> {
  typedef Departures<OBS>             Departures_;
  typedef DeparturesEnsemble<OBS>     DeparturesEnsemble_;
  typedef Geometry<MODEL>             Geometry_;
  typedef GeometryIterator<MODEL>     GeometryIterator_;
  typedef IncrementEnsemble4D<MODEL>  IncrementEnsemble4D_;
  typedef ObsErrors<OBS>              ObsErrors_;
  typedef ObsSpaces<OBS>              ObsSpaces_;

 public:
  static const std::string classname() {return "oops::LETKFSolverOOPS";}

  /// Constructor (allocates Wa, wa, saves options from the config)
  LETKFSolverOOPS(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t);

 private:
  /// Computes weights
  void computeWeights(const Departures_ &, const DeparturesEnsemble_ &,
                      const ObsErrors_ &) override;

  /// Applies weights and adds posterior inflation
  void applyWeights(const IncrementEnsemble4D_ &, IncrementEnsemble4D_ &,
                    const GeometryIterator_ &) override;

  LETKFSolverParameters options_;
  Eigen::MatrixXd Wa_;     // transformation matrix for ens. perts. Xa_=Xf*Wa
  Eigen::MatrixXd trans_;  // transformation matrix for ens. perts. Xa_=Xf*Wa+Xf*wa
  Eigen::VectorXd wa_;     // transformation matrix for ens. mean xa_=xf*wa
  Eigen::MatrixXd work_;   // temp array

  // eigen solver matrices
  Eigen::VectorXd eival_;
  Eigen::MatrixXd eivec_;

  const size_t nens_;      // ensemble size
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
LETKFSolverOOPS<MODEL, OBS>::LETKFSolverOOPS(ObsSpaces_ & obspaces, const Geometry_ & geom,
                                             const eckit::Configuration & config,
                                             size_t nens)
  : LETKFSolver<MODEL, OBS>(obspaces, geom, config, nens), nens_(nens)
{
  // read and parse options
  options_.deserialize(config);

  Log::info() << "Using OOPS implementation of LETKF" << std::endl;

  Log::info() << "Multiplicative inflation multCoeff=" <<
                    options_.infl.value().mult << std::endl;

  const LETKFInflationParameters & inflopt = options_.infl;
  if (inflopt.dortpp()) {
      Log::info() << "RTPP inflation will be applied with rtppCoeff=" <<
                      inflopt.rtpp << std::endl;
  } else {
      Log::info() << "RTPP inflation is not applied rtppCoeff is out of bounds (0,1],"
                  << " rtppCoeff=" << inflopt.rtpp << std::endl;
  }

  if (inflopt.dortps()) {
    Log::info() << "RTPS not implemented in LETKFOOPS solver"  << std::endl;
  }

  // pre-allocate transformation matrices
  trans_.resize(nens_, nens_);
  Wa_.resize(nens_, nens_);
  wa_.resize(nens_);
  work_.resize(nens_, nens_);

  // pre-allocate eigen sovler matrices
  eival_.resize(nens_);
  eivec_.resize(nens_, nens_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LETKFSolverOOPS<MODEL, OBS>::computeWeights(const Departures_ & dy,
                                            const DeparturesEnsemble_ & Yb,
                                            const ObsErrors_ & R) {
  // compute transformation matrix, save in Wa_, wa_, and trans_
  util::Timer timer(classname(), "computeWeights");

  const LETKFInflationParameters & inflopt = options_.infl;
  double infl = inflopt.mult;
  // fill in the work matrix (note that since the matrix is symmetric,
  // only lower triangular half is filled)
  // work = Y^T R^-1 Y + (nens_-1)/infl I
  for (size_t jj = 0; jj < nens_; ++jj) {
    Departures_ Cj(Yb[jj]);
    R.inverseMultiply(Cj);
    for (size_t ii = jj; ii < nens_; ++ii) {
      work_(ii, jj) = Cj.dot_product_with(Yb[ii]);
    }
    work_(jj, jj) += (nens_-1.0) / infl;
  }

  // eigenvalues and eigenvectors of the above matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(work_);
  eival_ = es.eigenvalues().real();
  eivec_ = es.eigenvectors().real();

  // Pa   = [ Yb^T R^-1 Yb + (nens_-1)/infl I ] ^-1
  work_ = eivec_ * eival_.cwiseInverse().asDiagonal() * eivec_.transpose();

  // Wa = sqrt[ (nens_-1) Pa ]
  trans_ = eivec_
    * ((nens_-1) * eival_.array().inverse()).sqrt().matrix().asDiagonal()
    * eivec_.transpose();

  // RTPP: Relaxation to prior perturbation.
  // RTPP is done on Wa before the ensemble mean translation is introduced
  if (inflopt.dortpp()) {
    trans_ = (1.-inflopt.rtpp)*trans_;
    trans_.diagonal() = Eigen::VectorXd::Constant(nens_, inflopt.rtpp)
                        + trans_.diagonal();
  }

  // wa = Pa Yb^T R^-1 dy
  Departures_ Rinvdy(dy);
  R.inverseMultiply(Rinvdy);
  for (size_t jj = 0; jj < nens_; ++jj) {
    wa_(jj) = Yb[jj].dot_product_with(Rinvdy);
  }
  wa_ = work_ * wa_;

  // add wa to each column of Wa to get T
  trans_.colwise() += wa_;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LETKFSolverOOPS<MODEL, OBS>::applyWeights(const IncrementEnsemble4D_ & bkg_pert,
                                          IncrementEnsemble4D_ & ana_pert,
                                          const GeometryIterator_ & i) {
// apply Wa_, wa_ (here combined in the trans_ matrix)
  util::Timer timer(classname(), "applyWeights");

  for (size_t itime=0; itime < bkg_pert[0].size(); ++itime) {
    for (size_t jj=0; jj < nens_; ++jj) {
      LocalIncrement gp = bkg_pert[0][itime].getLocal(i);
      gp *= trans_(0, jj);
      for (size_t ii=1; ii < nens_; ++ii) {
        LocalIncrement gp2 = bkg_pert[ii][itime].getLocal(i);
        gp2 *= trans_(ii, jj);
        gp += gp2;
      }
      ana_pert[jj][itime].setLocal(gp, i);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LETKFSOLVEROOPS_H_
