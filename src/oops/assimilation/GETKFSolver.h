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
#include "oops/assimilation/LETKFSolverParameters.h"
#include "oops/assimilation/LocalEnsembleSolver.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/IncrementEnsemble4D.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/State4D.h"
#include "oops/base/StateEnsemble4D.h"
#include "oops/generic/VerticalLocEV.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {

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
  typedef ObsDataVector<OBS, int>     ObsDataVector_;
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
  GETKFSolver(ObsSpaces_ &, const Geometry_ &, const eckit::Configuration &, size_t);

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
  LETKFSolverParameters options_;
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
                                const eckit::Configuration & config, size_t nens)
  : LocalEnsembleSolver<MODEL, OBS>(obspaces, geometry, config, nens),
    nens_(nens), geometry_(geometry),
    vertloc_(geometry_, config.getSubConfiguration("local ensemble DA.vertical localization")),
    neig_(vertloc_.neig()), nanal_(neig_*nens_), HZb_(obspaces, nanal_)
{
  options_.deserialize(config);
  const LETKFInflationParameters & inflopt = options_.infl;

  // parse inflation options
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
    Log::debug() << "Read H(X) from disk" << std::endl;
    // read modulated ensemble
    Observations_ ytmp(yb_mean);
    size_t ii = 0;
    for (size_t iens = 0; iens < nens_; ++iens) {
      for (size_t ieig = 0; ieig < neig_; ++ieig) {
        ytmp.read("hofxm"+std::to_string(iteration)+"_"+std::to_string(ieig+1)+
                      "_"+std::to_string(iens+1));
        HZb_[ii] = ytmp - yb_mean;
        Log::test() << "H(Zx) - ymean for member " << iens+1 << " eig "<< ieig+1
                    << " :" << std::endl << HZb_[ii] << std::endl;
        ii = ii + 1;
      }
    }
  } else {
    Log::debug() << "Computing H(X) online" << std::endl;
    // modulate ensemble of obs
    State4D_ xx_mean(ens_xx.mean());
    IncrementEnsemble4D_ dx(ens_xx, xx_mean, xx_mean[0].variables());
    IncrementEnsemble4D_ Ztmp(geometry_, xx_mean[0].variables(), ens_xx[0].validTimes(), neig_);
    size_t ii = 0;
    for (size_t iens = 0; iens < nens_; ++iens) {
      vertloc_.modulateIncrement(dx[iens], Ztmp);
      for (size_t ieig = 0; ieig < neig_; ++ieig) {
        State4D_ tmpState = xx_mean;
        tmpState += Ztmp[ieig];
        Observations_ tmpObs(this->obspaces_);
        this->computeHofX4D(tmpState, tmpObs);
        HZb_[ii] = tmpObs - yb_mean;
        tmpObs.save("hofxm"+std::to_string(iteration)+"_"+std::to_string(ieig+1)+
                      "_"+std::to_string(iens+1));
        Log::test() << "H(Zx) - ymean for member " << iens+1 << " eig "<< ieig+1 << " :"
                    << std::endl << HZb_[ii] << std::endl;
        ii = ii + 1;
      }
    }
  }
  this->readQC("EffectiveQC");
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
                 getkf_inflation, denkf, getkf);
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

  const LETKFInflationParameters & inflopt = options_.infl;

  // allocate tmp arrays
  LocalIncrement gptmpl = bkg_pert[0][0].getLocal(i);
  std::vector<double> tmp = gptmpl.getVals();
  size_t ngp = tmp.size();
  Eigen::MatrixXd Xb(ngp, nens_);

  // loop through analysis times and ens. members
  for (unsigned itime=0; itime < bkg_pert[0].size(); ++itime) {
    // cast bkg_pert ensemble at grid point i as an Eigen matrix Xb
    // modulates Xb
    Xb = vertloc_.modulateIncrement(bkg_pert, i, itime);

    // postmulptiply
    Eigen::VectorXd xa = Xb*wa_;  // ensemble mean update
    Eigen::MatrixXd Xa = Xb*Wa_;  // ensemble perturbation update

    // compute non-modulated Xb for RTPP and RTPS
    if (inflopt.dortpp() || inflopt.dortps()) {
      Xb.resize(ngp, nens_);
      for (size_t iens=0; iens < nens_; ++iens) {
        LocalIncrement gp = bkg_pert[iens][itime].getLocal(i);
        std::vector<double> tmp1 = gp.getVals();
        for (size_t iv=0; iv < ngp; ++iv) {
          Xb(iv, iens) = tmp1[iv];
        }
      }
    }

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

    // assign Xa_ to ana_pert
    for (size_t iens=0; iens < nens_; ++iens) {
      for (size_t iv=0; iv < ngp; ++iv) {
        tmp[iv] = Xa(iv, iens)+xa(iv);
      }
      gptmpl.setVals(tmp);
      ana_pert[iens][itime].setLocal(gptmpl, i);
    }
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
  std::vector<std::shared_ptr<ObsDataVector_>> outside;
  for (size_t jj = 0; jj < this->obspaces_.size(); ++jj) {
    outside.push_back(std::make_shared<ObsDataVector_>(this->obspaces_[jj],
            this->obspaces_[jj].obsvariables()));
  }
  locvector.ones();
  this->obsloc().computeLocalization(i, outside, locvector);
  locvector.mask(this->qcflags());
  Eigen::VectorXd local_omb_vec = this->omb_.packEigen(outside);

  if (local_omb_vec.size() == 0) {
    // no obs. so no need to update Wa_ and wa_
    // ana_pert[i]=bkg_pert[i]
    this->copyLocalIncrement(bkg_pert, i, ana_pert);
  } else {
    // if obs are present do normal KF update
    // get local Yb & HZ
    Eigen::MatrixXd local_Yb_mat = this->Yb_.packEigen(outside);
    Eigen::MatrixXd local_HZ_mat = this->HZb_.packEigen(outside);
    // create local obs errors
    Eigen::VectorXd local_invVarR_vec = this->invVarR_->packEigen(outside);
    // and apply localization
    Eigen::VectorXd localization = locvector.packEigen(outside);
    local_invVarR_vec.array() *= localization.array();
    computeWeights(local_omb_vec, local_HZ_mat, local_Yb_mat, local_invVarR_vec);
    applyWeights(bkg_pert, ana_pert, i);
  }
}

}  // namespace oops
#endif  // OOPS_ASSIMILATION_GETKFSOLVER_H_
