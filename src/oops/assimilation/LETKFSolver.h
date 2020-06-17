/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

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


#ifndef OOPS_ASSIMILATION_LETKFSOLVER_H_
#define OOPS_ASSIMILATION_LETKFSOLVER_H_

#include <Eigen/Dense>
#include <cfloat>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/LETKFSolverBase.h"
#include "oops/assimilation/LETKFSolverParameters.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/ObsErrors.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL>
class LETKFSolver : public LETKFSolverBase<MODEL> {
  typedef Departures<MODEL>         Departures_;
  typedef DeparturesEnsemble<MODEL> DeparturesEnsemble_;
  typedef GeometryIterator<MODEL>   GeometryIterator_;
  typedef IncrementEnsemble<MODEL>  IncrementEnsemble_;
  typedef ObsErrors<MODEL>          ObsErrors_;

 public:
  LETKFSolver(const eckit::Configuration &, size_t);

 private:
  /// Computes weights
  void computeWeights(const Departures_ &, const DeparturesEnsemble_ &,
                      const ObsErrors_ &) override;

  /// Applies weights and adds posterior inflation
  void applyWeights(const IncrementEnsemble_ &, IncrementEnsemble_ &,
                    const GeometryIterator_ &) override;

  LETKFSolverParameters options_;
  Eigen::MatrixXd Wa_;  // transformation matrix for ens. perts. Xa=Xf*Wa
  Eigen::VectorXd wa_;  // transformation matrix for ens. mean xa=xf*wa

  // eigen solver matrices
  Eigen::VectorXd eival_;
  Eigen::MatrixXd eivec_;

  // parameters
  int nens_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LETKFSolver<MODEL>::LETKFSolver(const eckit::Configuration & config,
                                size_t nens)
                                : nens_(nens) {
  options_.deserialize(config);
  const LETKFInflationParameters & inflopt = options_.infl;

  Log::info() << "Using EIGEN implementation of LETKF" << std::endl;

  Log::info() << "Multiplicative inflation multCoeff=" <<
                 inflopt.mult << std::endl;

  if (inflopt.dortpp()) {
      Log::info() << "RTPP inflation will be applied with rtppCoeff=" <<
                    inflopt.rtpp << std::endl;
  } else {
      Log::info() << "RTPP inflation is not applied rtppCoeff is out of bounds (0,1], rtppCoeff="
                  << inflopt.rtpp << std::endl;
  }

  if (inflopt.dortps()) {
    Log::info() << "RTPS inflation will be applied with rtpsCoeff=" <<
                    inflopt.rtps << std::endl;
  } else {
    Log::info() << "RTPS inflation is not applied rtpsCoeff is out of bounds (0,1], rtpsCoeff="
                << inflopt.rtps << std::endl;
  }

  // pre-allocate transformation matrices
  Wa_.resize(nens_, nens_);
  wa_.resize(nens_);

  // pre-allocate eigen sovler matrices
  eival_.resize(nens_);
  eivec_.resize(nens_, nens_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LETKFSolver<MODEL>::computeWeights(const Departures_ & dy_oops,
                                        const DeparturesEnsemble_ & Yb_oops,
                                        const ObsErrors_ & R_oops) {
  // compute transformation matrix, save in Wa_, wa_
  // uses C++ eigen interface
  // implements LETKF from Hunt et al. 2007
  const LETKFInflationParameters & inflopt = options_.infl;

  // cast oops objects to eigen
  Eigen::MatrixXd dy = dy_oops.packEigen();

  Eigen::MatrixXd Yb = Yb_oops.packEigen();
  Eigen::MatrixXd diagInvR = R_oops.packInverseVarianceEigen();

  // fill in the work matrix
  // work = Y^T R^-1 Y + (nens_-1)/infl I
  double infl = inflopt.mult;
  Eigen::MatrixXd work = Yb*(diagInvR.asDiagonal()*Yb.transpose());
  work.diagonal() += Eigen::VectorXd::Constant(nens_, (nens_-1)/infl);

  // eigenvalues and eigenvectors of the above matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(work);
  eival_ = es.eigenvalues().real();
  eivec_ = es.eigenvectors().real();

  // Pa   = [ Yb^T R^-1 Yb + (nens_-1)/infl I ] ^-1
  work = eivec_ * eival_.cwiseInverse().asDiagonal() * eivec_.transpose();

  // Wa = sqrt[ (nens_-1) Pa ]
  Wa_ = eivec_
      * ((nens_-1) * eival_.array().inverse()).sqrt().matrix().asDiagonal()
      * eivec_.transpose();

  // wa = Pa Yb^T R^-1 dy
  wa_ = work * (Yb * (diagInvR.asDiagonal()*dy.transpose()));
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LETKFSolver<MODEL>::applyWeights(const IncrementEnsemble_ & bkg_pert,
                                      IncrementEnsemble_ & ana_pert,
                                      const GeometryIterator_ & i) {
  // applies Wa_, wa_

  const LETKFInflationParameters & inflopt = options_.infl;

  LocalIncrement gptmpl = bkg_pert[0][0].getLocal(i);
  std::vector<double> tmp1 = gptmpl.getVals();
  size_t ngp = tmp1.size();

  // loop through analysis times and ens. members
  for (size_t itime=0; itime < bkg_pert[0].size(); ++itime) {
    // make grid point forecast pert ensemble array
    Eigen::MatrixXd Xb(ngp, nens_);
    for (int iens=0; iens < nens_; ++iens) {
      LocalIncrement gp = bkg_pert[iens][itime].getLocal(i);
      std::vector<double> tmp = gp.getVals();
      for (size_t iv=0; iv < ngp; ++iv) {
        Xb(iv, iens) = tmp[iv];
      }
    }

    // postmulptiply
    Eigen::VectorXd xa = Xb*wa_;   // ensemble mean update
    Eigen::MatrixXd Xa = Xb*Wa_;   // ensemble perturbation update

    // RTPP inflation
    if (inflopt.dortpp()) {
      Xa = (1-inflopt.rtpp)*Xa+inflopt.rtpp*Xb;
    }

    // RTPS inflation
    double eps = DBL_EPSILON;
    if (inflopt.dortps()) {
      // posterior spread
      Eigen::ArrayXd asprd = Xa.array().square().rowwise().sum()/(nens_-1);
      asprd = asprd.sqrt();
      asprd = (asprd < eps).select(eps, asprd);  // avoid nan overflow for vars with no spread

      // prior spread
      Eigen::ArrayXd fsprd = Xb.array().square().rowwise().sum()/(nens_-1);
      fsprd = fsprd.sqrt();
      fsprd = (fsprd < eps).select(eps, fsprd);

      // rtps inflation factor
      Eigen::ArrayXd rtpsInfl = inflopt.rtps*((fsprd-asprd)/asprd) + 1;
      rtpsInfl = (rtpsInfl < inflopt.rtpsInflMin()).select(inflopt.rtpsInflMin(), rtpsInfl);
      rtpsInfl = (rtpsInfl > inflopt.rtpsInflMax()).select(inflopt.rtpsInflMax(), rtpsInfl);

      // inlfate perturbation matrix
      Xa.array().colwise() *= rtpsInfl;
    }

    // assign Xa to ana_pert
    for (int iens=0; iens < nens_; ++iens) {
      for (size_t iv=0; iv < ngp; ++iv) {
        tmp1[iv] = Xa(iv, iens)+xa(iv);   // if Xa = Xb*Wa;
      }
      gptmpl.setVals(tmp1);
      ana_pert[iens][itime].setLocal(gptmpl, i);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LETKFSOLVER_H_
