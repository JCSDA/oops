/*
 * (C) Copyright 2020-2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_LOCALENSEMBLESOLVER_H_
#define OOPS_ASSIMILATION_LOCALENSEMBLESOLVER_H_

#include <Eigen/Dense>
#include <algorithm>
#include <cfloat>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/assimilation/DFSCalculator.h"
#include "oops/assimilation/SubensembleSplitter.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/Geometry.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/LinearModel.h"
#include "oops/base/Model.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObserversTLAD.h"
#include "oops/base/ObsLocalizations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/generic/PseudoLinearModelIncrement4D.h"
#include "oops/generic/PseudoModelState4D.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace oops {
  class Variables;

/// \brief Base class for local ensemble solvers
template <typename MODEL, typename OBS>
class LocalEnsembleSolver {
  typedef Observers<MODEL, OBS>       Observers_;
  typedef ObserversTLAD<MODEL, OBS>   ObserversTLAD_;
  typedef Departures<OBS>             Departures_;
  typedef DeparturesEnsemble<OBS>     DeparturesEnsemble_;
  typedef Geometry<MODEL>             Geometry_;
  typedef GeometryIterator<MODEL>     GeometryIterator_;
  typedef IncrementSet<MODEL>         IncrementSet_;
  typedef ObsAuxControls<OBS>         ObsAux_;
  typedef ObsAuxIncrements<OBS>       ObsAuxInc_;
  typedef ObsDataVector<OBS, int>     ObsDataInt_;
  typedef ObsEnsemble<OBS>            ObsEnsemble_;
  typedef ObsError<OBS>               ObsError_;
  typedef ObsErrors<OBS>              ObsErrors_;
  typedef Observations<OBS>           Observations_;
  typedef ObsLocalizations<MODEL, OBS> ObsLocalizations_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef PseudoModelState4D<MODEL>   PseudoModel_;
  typedef PseudoLinearModelIncrement4D<MODEL> PseudoLinearModel_;
  typedef State<MODEL>                State_;
  typedef StateSet<MODEL>             StateSet_;
  typedef State4D<MODEL>              State4D_;
  typedef Increment<MODEL>            Increment_;
  typedef Increment4D<MODEL>          Increment4D_;
  typedef LinearModel<MODEL>          LinearModel_;
  typedef Model<MODEL>                Model_;
  typedef ModelAuxControl<MODEL>      ModelAux_;
  typedef ModelAuxIncrement<MODEL>    ModelAuxInc_;
  typedef ObsDataVector<OBS, int>     ObsData_;
  typedef std::vector<std::shared_ptr<ObsData_>> ObsDataVec_;

 public:
  static const std::string classname() {return "oops::LocalEnsembleSolver";}

  /// initialize solver with \p obspaces, \p geometry, full \p config and \p nens ensemble size
  /// \p xbmean state is used if an implementation needs a reference state
  /// solver will use a list of analysis variables specified in \p incvars
  LocalEnsembleSolver(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                      const eckit::Configuration & config, size_t nens, const StateSet_ & xbmean,
                      const Variables & incvars);
  virtual ~LocalEnsembleSolver() = default;

  /// update background ensemble \p bg to analysis ensemble \p an for all points on this PE
  virtual void measurementUpdate(const IncrementSet_ & bg,
                                 IncrementSet_ & an);

  /// update background ensemble \p bg to analysis ensemble \p an at a grid point location \p i
  virtual void measurementUpdate(const Eigen::VectorXd &,
                                 const ObsErrors_ &,
                                 const Departures_ &,
                                 const IncrementSet_ &,
                                 const GeometryIterator_ &,
                                 IncrementSet_ &) = 0;

  /// copy \p an[\p i] = \p bg[\p i] (e.g. when there are no local observations to update state)
  virtual void copyLocalIncrement(const IncrementSet_ & bg,
                                  const GeometryIterator_ & i,
                                  IncrementSet_ & an) const;

  /// apply posterior inflation to a local ensemble
  void posteriorInflation(const Eigen::MatrixXd & Xb, Eigen::MatrixXd & Xa) const;

  /// accessor to obs localizations
  const ObsLocalizations_ & obsloc() const {return *obsloc_;}

  /// Enable Nerger et al. 2012 observation localization regulation
  bool useNergerRegulation() const { return useNergerRegulation_; }

 protected:
  const Geometry_  & geometry_;   ///< Geometry associated with the updated states
  const ObsSpaces_ & obspaces_;   ///< ObsSpaces used in the update
  Observations_ ybmean_;
  Departures_ omb_;               ///< obs - mean(H(x)); set in computeHofX method
  std::unique_ptr<DeparturesEnsemble_> Yb_;   ///< ensemble perturbations in the observation space;
                                              ///< set in computeHofX method
  std::unique_ptr<ObsErrors_>  R_;         ///< observation errors, set in computeHofX method
  std::unique_ptr<Departures_> invVarR_;   ///< inverse observation error variance for assimilated
                                           ///< observations
  eckit::LocalConfiguration inflopt_;
  int unperturbedIdx_;  ///< For stochastic filters: ensemble member index where obs are
                        ///< (optionally) unperturbed

  const StateSet_ & xbmean_;     ///< ensemble mean or a control member that will be used to
                                 ///< center the prior ensemble
  const Variables incvars_;
  const eckit::LocalConfiguration obsconf_;  ///< configuration for observations
  const eckit::LocalConfiguration observersconf_;  ///< configuration for observations.observers
  std::unique_ptr<ObsLocalizations_> obsloc_;          ///< observation space localization

  bool doCrossValidation;  ///< cross validation trigger
  size_t nsubens_;  ///< no. of subensembles
  std::unique_ptr<oops::SubensembleSplitter> SubensembleSplitter_;  ///< pointer to splitter

  std::unique_ptr<DFSCalculator<OBS>> dfsCalculator_;  /// DFS
  bool useNergerRegulation_;  ///< toggle for Nerger observation localisation regulation

  void applyAssimilatedMask(Departures_ & dep) const { dep.mask(*invVarR_); }
  /// Return true if the config requests SVD for the inverse analysis error covariance
  /// decomposition, false if eigendecomposition is requested. The string can be either
  /// "singular value decomposition" or "svd" (case insensitive) for SVD, and "eigendecomposition"
  /// for eigendecomposition (also case insensitive and default).
  const bool svdRequested(const eckit::Configuration & config) {
    std::string method =
        config.getString("local ensemble DA.inverse analysis error covariance decomposition method",
                         "eigendecomposition");
    std::transform(method.begin(), method.end(), method.begin(), ::tolower);
    // Assert all alphabetic characters in string are lowercase
    assert(std::all_of(method.begin(), method.end(),
                       [](unsigned char c) { return !std::isalpha(c) || std::islower(c); }));
    return (method == "singular value decomposition" || method == "svd");
  }

 private:
  // compute local inverse R vector with optional Nerger regulation
  Departures_ computeNergerLocalR(const Departures_ & locvector) const;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
LocalEnsembleSolver<MODEL, OBS>::LocalEnsembleSolver(ObsSpaces_ & obspaces,
                                                     const Geometry_ & geometry,
                                                     const eckit::Configuration & config,
                                                     size_t nens,
                                                     const StateSet_ & xbmean,
                                                     const Variables & incvars)
  : geometry_(geometry),
    obspaces_(obspaces),
    omb_(obspaces_),
    ybmean_(obspaces_),
    inflopt_(config.getSubConfiguration("local ensemble DA.inflation")),
    xbmean_(xbmean),
    incvars_(incvars),
    obsconf_(config, "observations"),
    observersconf_(obsconf_, "observers"),
    nsubens_(1) {
  // initialize and print options

  if (config.has("local ensemble DA.unperturbed obs ensemble member index")) {
    unperturbedIdx_ = config.getInt("local ensemble DA.unperturbed obs ensemble member index");
    if (unperturbedIdx_ < 0 || unperturbedIdx_ >= static_cast<int>(nens)) {
      const std::string e = "unperturbed obs ensemble member index out of bounds " +
                            std::to_string(unperturbedIdx_) + " for ensemble size " +
                            std::to_string(nens) + "\nMust be in the closed interval [0, " +
                            std::to_string(nens - 1) + "]";
      oops::Log::error() << e << std::endl;
      throw eckit::BadParameter(e);
    }
  } else {
    unperturbedIdx_ = -1;  // no unperturbed member
  }
  useNergerRegulation_ = config.getBool("local ensemble DA.use nerger regulation", false);
  Log::info() << "Nerger et al. (2012) observation localization regulation: " <<
                 (useNergerRegulation_ ? "ON" : "OFF (default)") << std::endl;
  Log::info() << "Multiplicative inflation will be applied with multCoeff=" <<
                 inflopt_.getDouble("mult", 1.0) << std::endl;
  const double rtpp = inflopt_.getDouble("rtpp", 0.0);
  const double rtps = inflopt_.getDouble("rtps", 0.0);
  if (rtpp > 0.0 && rtpp <= 1.0) {
      Log::info() << "RTPP inflation will be applied with rtppCoeff=" <<
                    rtpp << std::endl;
  } else {
      Log::info() << "RTPP not applied: rtppCoeff is out of bounds (0,1], rtppCoeff="
                  << rtpp << std::endl;
  }
  if (rtps > 0.0) {
    Log::info() << "RTPS inflation will be applied with rtpsCoeff=" <<
                    rtps << std::endl;
  } else {
    Log::info() << "RTPS not applied: rtpsCoeff is <=0, rtpsCoeff="
                << rtps << std::endl;
  }

  // Instantiate the DFS calculator
  if (config.has("driver")) {
    const eckit::LocalConfiguration driver(config, "driver");
    if (driver.getBool("dfs", false)) {
      const size_t commSize = geometry_.getComm().size();
      if (commSize != 1) {
        const std::string e =
          "DFS estimator is not yet implemented for observations distributed over multiple MPI "
          "tasks.  Run with one task or set 'driver.dfs' to false.";
        oops::Log::error() << e << std::endl;
        throw eckit::NotImplemented(e, Here());
      }
      dfsCalculator_ = std::make_unique<DFSCalculator<OBS>>(omb_.nobs());
    }
  }

  ObsEnsemble_ obsens(obspaces_, nens);
  Observations_ y_mean_xb(obspaces_);
  for (size_t jj = 0; jj < nens; ++jj) {
    obsens[jj].read("hofx0_"+std::to_string(jj+1));
  }
  y_mean_xb.read("hofx_y_mean_xb0");
  ybmean_ = obsens.mean();
  Departures_ tmpDeps(this->obspaces_);
  // initialize Yb_ and obsloc_
  Yb_ = std::make_unique<DeparturesEnsemble_>(obspaces_, nens);
  obsloc_ = std::make_unique<ObsLocalizations_>(observersconf_, obspaces_);
  for (size_t iens = 0; iens < Yb_->size(); ++iens) {
    tmpDeps = obsens[iens] - ybmean_;
    // updateAssimilatedMask(tmpDeps);
    Yb_->setData(iens, tmpDeps);
  }
  Observations_ yobs(obspaces_, "ObsValue");
  omb_ = yobs - ybmean_;
  R_.reset(new ObsErrors_(this->observersconf_, this->obspaces_));
  Observations_ obserr(obspaces_, "EffectiveError0");
  for (size_t jj = 0; jj < R_->size(); ++jj) {
    (*R_)[jj].update(obserr[jj]);
  }
  invVarR_ = std::make_unique<Departures_>(R_->inverseVariance());
  oops::Log::trace() << "LocalEnsembleSolver created" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
typename LocalEnsembleSolver<MODEL, OBS>::Departures_
LocalEnsembleSolver<MODEL, OBS>::computeNergerLocalR(
  const Departures_ & locvector) const {
  // Nerger et al. (2012) observation localization regulation applied to R matrix
  // Currently only supports diagonal R

  Log::trace() << "LocalEnsembleSolver<MODEL, OBS>::computeNergerLocalR starting" << std::endl;

  if (!invVarR_) {
    oops::Log::error() << "invVarR_ not initialized before computeNergerLocalR. \n"
                       << "Fail because this implies that observation error is a \n"
                       << "correlated matrix not diagonal" << std::endl;
    throw std::logic_error("invVarR_ is null");
  }
  const Eigen::VectorXd invVarR_local = invVarR_->packEigen(locvector);
  const Eigen::MatrixXd Yb_local = this->Yb_->packEigen(locvector).template cast<double>();

  // VarP_i = (1/(nens-1)) * sum_m Yb[m]_i^2, or the diagonal of background covariance error
  // per local observation
  const Eigen::ArrayXd VarP_local = Yb_local.array().square().colwise().sum().transpose() /
                                     static_cast<double>(Yb_local.rows() - 1);

  // gammaNerger_i = VarP_i / VarR_i = VarP_i * invVarR_i, per local observation
  const Eigen::ArrayXd gammaNerger_local = VarP_local * invVarR_local.array();

  // betaNerger_i = 1 + gammaNerger_i * (1 - loc_i), per local observation
  const Eigen::ArrayXd loc_local = locvector.packEigen(locvector).array();
  const Eigen::ArrayXd betaNerger_arr = 1.0 + gammaNerger_local * (1.0 - loc_local);

  // Build nergerBeta Departures_: 1 at non-local obs, betaNerger_i at local obs positions
  Departures_ nergerBeta(obspaces_);
  nergerBeta.ones();
  const std::vector<std::vector<size_t>> localIndices =
      nergerBeta.maskAndSerialIndices(locvector);
  const std::vector<size_t> serialSizes = nergerBeta.serialSizes();
  std::vector<double> betaVec;
  betaVec.reserve(nergerBeta.serialSize());
  for (size_t ii = 0; ii < nergerBeta.size(); ++ii) {
    nergerBeta[ii].serialize(betaVec);
  }
  size_t iStart = 0;
  size_t iObs = 0;
  for (size_t ii = 0; ii < localIndices.size(); ++ii) {
    for (const size_t idx : localIndices[ii]) {
      betaVec[idx + iStart] = betaNerger_arr(iObs++);
    }
    iStart += serialSizes[ii];
  }
  ASSERT(iObs == static_cast<size_t>(betaNerger_arr.size()));
  size_t deserIdx = 0;
  for (size_t ii = 0; ii < nergerBeta.size(); ++ii) {
    nergerBeta[ii].deserialize(betaVec, deserIdx);
  }
  return nergerBeta;
}

template <typename MODEL, typename OBS>
void LocalEnsembleSolver<MODEL, OBS>::measurementUpdate
(const IncrementSet_ & bkg_pert, IncrementSet_ & ana_pert) {
  for (GeometryIterator_ i = geometry_.begin(); i != geometry_.end(); ++i) {
    // create the local subset of observations
    Departures_ locvector(this->obspaces_);
    locvector.ones();
    this->obsloc().computeLocalization(i, locvector);
    this->applyAssimilatedMask(locvector);
    const Eigen::VectorXd local_omb_vec = this->omb_.packEigen(locvector);
    if (local_omb_vec.size() == 0) {
      // no obs. so no need to update Wa_ and wa_
      // ana_pert[i] = bkg_pert[i]
      this->copyLocalIncrement(bkg_pert,
                               i,
                               ana_pert);
    } else {
      if (useNergerRegulation_) {
        const Departures_ localizationMask(locvector);
        const Departures_ nergerBeta = this->computeNergerLocalR(locvector);
        ASSERT(locvector.serialSize() == nergerBeta.serialSize());
        locvector /= nergerBeta;
        locvector.mask(localizationMask);
    }
    R_->localize(locvector);

    this->measurementUpdate(local_omb_vec,
                            *R_,
                            locvector,
                            bkg_pert,
                            i,
                            ana_pert);
    }
  }
  // Calculate the DFS >>>
  if (dfsCalculator_) {
    dfsCalculator_->finalize();
    Departures_ mask(this->obspaces_);
    mask.ones();
    this->applyAssimilatedMask(mask);
    // DFS per observation block (save in the log file and print on the screen)
    const auto stats = dfsCalculator_->computeBlockStats(omb_, mask, observersconf_);
    dfsCalculator_->printBlockStats(stats);
    // DFS per observation (save in nc file)
    dfsCalculator_->savePerObservation(omb_, mask, this->obspaces_, "DFS");
  }
}
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LocalEnsembleSolver<MODEL, OBS>::copyLocalIncrement(const IncrementSet_ & bkg_pert,
                                                         const GeometryIterator_ & i,
                                                         IncrementSet_ & ana_pert) const {
  // ana_pert[i]=bkg_pert[i]
  for (size_t itime=0; itime < bkg_pert.time_size(); ++itime) {
    for (size_t iens=0; iens < bkg_pert.ens_size(); ++iens) {
      LocalIncrement gp = bkg_pert(itime, iens).getLocal(i);
      ana_pert(itime, iens).setLocal(gp, i);
    }
  }
}

// -----------------------------------------------------------------------------

/// \brief factory for LocalEnsembleSolver solvers
template <typename MODEL, typename OBS>
class LocalEnsembleSolverFactory {
  typedef Geometry<MODEL>           Geometry_;
  typedef ObsSpaces<OBS>            ObsSpaces_;
  typedef StateSet<MODEL>           StateSet_;
 public:
  static std::unique_ptr<LocalEnsembleSolver<MODEL, OBS>> create(ObsSpaces_ &, const Geometry_ &,
                                                        const eckit::Configuration &,
                                                        size_t, const StateSet_ &,
                                                        const Variables &);
  virtual ~LocalEnsembleSolverFactory() = default;
 protected:
  explicit LocalEnsembleSolverFactory(const std::string &);
 private:
  virtual LocalEnsembleSolver<MODEL, OBS> * make(ObsSpaces_ &, const Geometry_ &,
                                        const eckit::Configuration &, size_t,
                                        const StateSet_ &, const Variables &) = 0;
  static std::map < std::string, LocalEnsembleSolverFactory<MODEL, OBS> * > & getMakers() {
    static std::map < std::string, LocalEnsembleSolverFactory<MODEL, OBS> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class OBS, class T>
class LocalEnsembleSolverMaker : public LocalEnsembleSolverFactory<MODEL, OBS> {
  typedef Geometry<MODEL>           Geometry_;
  typedef ObsSpaces<OBS>            ObsSpaces_;
  typedef StateSet<MODEL>           StateSet_;

  virtual LocalEnsembleSolver<MODEL, OBS> * make(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                        const eckit::Configuration & conf, size_t nens,
                                        const StateSet_ & xbmean, const Variables & incvars)
    { return new T(obspaces, geometry, conf, nens, xbmean, incvars); }
 public:
  explicit LocalEnsembleSolverMaker(const std::string & name)
    : LocalEnsembleSolverFactory<MODEL, OBS>(name) {}
};

// =============================================================================

template <typename MODEL, typename OBS>
LocalEnsembleSolverFactory<MODEL, OBS>::LocalEnsembleSolverFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in local ensemble solver factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<LocalEnsembleSolver<MODEL, OBS>>
LocalEnsembleSolverFactory<MODEL, OBS>::create(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                                  const eckit::Configuration & conf, size_t nens,
                                  const StateSet_ & xbmean, const Variables & incvars) {
  Log::trace() << "LocalEnsembleSolver<MODEL, OBS>::create starting" << std::endl;
  const std::string id = conf.getString("local ensemble DA.solver");
  typename std::map<std::string, LocalEnsembleSolverFactory<MODEL, OBS>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in local ensemble solver factory." << std::endl;
    Log::error() << "Local ensemble solver factory has "
                 << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, LocalEnsembleSolverFactory<MODEL, OBS>*>::const_iterator
         jj = getMakers().begin(); jj != getMakers().end(); ++jj) {
       Log::error() << "A " << jj->first << " LocalEnsembleSolver" << std::endl;
    }
    throw std::runtime_error(id + " does not exist in local ensemble solver factory.");
  }
  std::unique_ptr<LocalEnsembleSolver<MODEL, OBS>>
    ptr(jloc->second->make(obspaces, geometry, conf, nens, xbmean, incvars));
  Log::trace() << "LocalEnsembleSolver<MODEL, OBS>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LOCALENSEMBLESOLVER_H_
