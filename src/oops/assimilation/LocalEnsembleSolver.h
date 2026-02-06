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
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
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
#include "oops/base/StateEnsemble4D.h"
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
  typedef StateEnsemble4D<MODEL>      StateEnsemble4D_;
  typedef PseudoModelState4D<MODEL>   PseudoModel_;
  typedef PseudoLinearModelIncrement4D<MODEL> PseudoLinearModel_;
  typedef State<MODEL>                State_;
  typedef StateSet<MODEL>             StateSet_;
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

  /// computes ensemble H(\p xx), returns mean H(\p xx), saves as hofx \p iteration
  virtual Observations_ computeHofX(const StateEnsemble4D_ & xx, size_t iteration,
                      bool readFromDisk);

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
  bool useLinearObserver() const { return useLinearObserver_; }

 protected:
  const Geometry_  & geometry_;   ///< Geometry associated with the updated states
  const ObsSpaces_ & obspaces_;   ///< ObsSpaces used in the update
  Departures_ omb_;               ///< obs - mean(H(x)); set in computeHofX method
  std::unique_ptr<DeparturesEnsemble_> Yb_;   ///< ensemble perturbations in the observation space;
                                              ///< set in computeHofX method
  std::unique_ptr<ObsErrors_>  R_;         ///< observation errors, set in computeHofX method
  std::unique_ptr<Departures_> invVarR_;   ///< inverse observation error variance for assimilated
                                           ///< observations; set in initializeAssimilatedMask
  std::vector<ObsDataInt_> qcflags_;  ///< quality control flags
  std::unique_ptr<ObserversTLAD_> linear_hofx_;  ///< linear observer
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

  /// Create a mask that excludes observations which will not be assimilated (e.g. failed QC) using
  /// a single ensemble member.
  void initializeAssimilatedMask() {
    // Inverse variances have missing values where obs have failed QC, though R_
    // is only valid for a single ensemble member
    invVarR_ = std::make_unique<Departures_>(R_->inverseVariance());
  }
  /// Update departures for assimilated observations by masking with \p mask
  /// (e.g. nonzero QC flags or missing obs) - should only be done in the
  /// computeHofX method.
  void updateAssimilatedMask(const Departures_ & mask) { invVarR_->mask(mask); }
  /// Apply the assimilated mask to a departures vector \p dep - missing values
  /// in the mask will become missing in the departures.
  void applyAssimilatedMask(Departures_ & dep) const { dep.mask(*invVarR_); }
  /// Apply the non-linear observation operator to the background state \p xx and store the result
  /// in \p yy. If \ref useLinearObserver() returns true, the observation operator is also
  /// linearized about the background state \p xx (cf \ref applyLinearToPerturbations).
  void computeHofX4D(const eckit::Configuration & config, const StateSet_ & xx, Observations_ & yy,
                     const util::Duration & flength, const util::Duration & default_tstep,
                     const ObsAux_ & obsaux, const ModelAux_ & moderr,
                     ObsErrors_ & R, std::vector<ObsDataInt_> & qcflags);
  /// Runs a linear model on 4D perturbations from the ensemble mean ( \p xx - \ref xbmean_ ) and
  /// applies a linearized observation operator \ref linear_hofx_ to background departures in the
  /// ensemble \ref Yb_ at ensemble index \p iens .
  /// Intended for following the procedure from Shlyaeva, A., & Whitaker, J. S. (2018)
  /// https://doi.org/10.1029/2018MS001309 where the linearized observation operator is applied to
  /// ensemble perturbations. In more detail, the obs operator is linearized about the mean state
  /// (H^tilde in eq. 6) - this happens in \ref computeHofX4D when \ref useLinearObserver() returns
  /// true. This is applied to the ensemble perturbations from the mean (X) (Eq 6.). The result is
  /// added to the application of the nonlinear observation operator to the mean state (Eq. 5).
  void applyLinearToPerturbations(const Increment4D_ & dx, const util::Duration & flength,
                                  const util::Duration & default_tstep,
                                  const ObsAuxInc_ & obsauxinc, const ModelAuxInc_ & moderrinc,
                                  const PostProcessor<Increment_> & posttl,
                                  const PostProcessorTLAD<MODEL> & posttrajtl,
                                  Departures_ & tmpDeps);
  /// Read pre-calculated HofX \p y_mean_xb from disk for each observation ensemble member \p obsens
  /// of size \p nens at \p iteration.
  void readHofX(ObsEnsemble_ & obsens, const size_t nens, const size_t iteration,
                Observations_ & y_mean_xb);
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
  bool useLinearObserver_;
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
    inflopt_(config.getSubConfiguration("local ensemble DA.inflation")),
    xbmean_(xbmean),
    incvars_(incvars),
    obsconf_(config, "observations"),
    observersconf_(obsconf_, "observers"),
    nsubens_(1) {
  // initialize and print options

  useLinearObserver_ = config.getBool("local ensemble DA.use linear observer", false);
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
  Log::info() << "Multiplicative inflation will be applied with multCoeff=" <<
                 inflopt_.getDouble("mult", 1.0) << std::endl;
  const double rtpp = inflopt_.getDouble("rtpp", 0.0);
  const double rtps = inflopt_.getDouble("rtps", 0.0);
  if (rtpp > 0.0 && rtpp <= 1.0) {
      Log::info() << "RTPP inflation will be applied with rtppCoeff=" <<
                    rtpp << std::endl;
  } else {
      Log::info() << "RTPP inflation is not applied rtppCoeff is out of bounds (0,1], rtppCoeff="
                  << rtpp << std::endl;
  }
  if (rtps > 0.0 && rtps <= 1.0) {
    Log::info() << "RTPS inflation will be applied with rtpsCoeff=" <<
                    rtps << std::endl;
  } else {
    Log::info() << "RTPS inflation is not applied rtpsCoeff is out of bounds (0,1], rtpsCoeff="
                << rtps << std::endl;
  }
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    ObsDataInt_ qcflags(obspaces_[jj], obspaces_[jj].obsvariables());
    qcflags_.push_back(qcflags);
  }
}

// -----------------------------------------------------------------------------

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
    const Eigen::VectorXd localization = locvector.packEigen(locvector);
    R_->localize(locvector);

    if (local_omb_vec.size() == 0) {
      // no obs. so no need to update Wa_ and wa_
      // ana_pert[i] = bkg_pert[i]
      this->copyLocalIncrement(bkg_pert,
                               i,
                               ana_pert);
    } else {
      this->measurementUpdate(local_omb_vec,
                              *R_,
                              locvector,
                              bkg_pert,
                              i,
                              ana_pert);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LocalEnsembleSolver<MODEL, OBS>::computeHofX4D(const eckit::Configuration & config,
                                                    const StateSet_ & xx, Observations_ & yy,
                                                    const util::Duration & flength,
                                                    const util::Duration & default_tstep,
                                                    const ObsAux_ & obsaux,
                                                    const ModelAux_ & moderr,
                                                    ObsErrors_ & Rmat,
                                                    std::vector<ObsDataInt_> & qcflags) {
  // Setup pseudo model to run on ensemble mean
  State_ init_xx = xx[0];
  std::unique_ptr<PseudoModel_> pseudomodel(new PseudoModel_(xx, default_tstep));
  const Model_ model(std::move(pseudomodel));

  // setup nonlinear postprocessor nonlinear observers
  PostProcessor<State_> post;
  Observers_ hofx(this->obspaces_, this->obsconf_);

  // initialize nonlinear model postprocessor
  hofx.initialize(this->geometry_, obsaux, Rmat, post, config);

  if (useLinearObserver()) {
    // Set up linear observer
    linear_hofx_ = std::make_unique<ObserversTLAD_>(obspaces_, obsconf_);
    // Set up linear postprocessor
    PostProcessorTLAD<MODEL> posttraj;
    // add linearized H(x) to the nonlinear model postprocessor
    linear_hofx_->initializeTraj(this->geometry_, obsaux, posttraj);
    // create TrajectorySaver with hofx_linear, and enroll in post
    post.enrollProcessor(new TrajectorySaver<MODEL>(eckit::LocalConfiguration(),
                                                    this->geometry_, posttraj));
  }

  // run nonlinear model on the ensemble mean
  model.forecast(init_xx, moderr, flength, post);

  // compute nonlinear H(x)
  hofx.finalize(yy, qcflags);

  if (useLinearObserver()) {
    // set the background state to linearise H about.
    linear_hofx_->finalizeTraj(qcflags);
  }
}
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LocalEnsembleSolver<MODEL, OBS>::applyLinearToPerturbations(
    const Increment4D_ & dx, const util::Duration & flength, const util::Duration & default_tstep,
    const ObsAuxInc_ & obsauxinc, const ModelAuxInc_ & moderrinc,
    const PostProcessor<Increment_> & posttl, const PostProcessorTLAD<MODEL> & posttrajtl,
    Departures_ & tmpDeps) {
  Increment_ init_dx = dx[0];
  std::unique_ptr<PseudoLinearModel_> pseudolinearmodel =
        std::make_unique<PseudoLinearModel_>(dx, default_tstep);
  const LinearModel_ linear_model(std::move(pseudolinearmodel));
  // run linear model on the ensemble perturbation, compute linear H*dx
  linear_model.forecastTL(init_dx, moderrinc, flength, posttl, posttrajtl);
  linear_hofx_->finalizeTL(obsauxinc, tmpDeps);
}
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LocalEnsembleSolver<MODEL, OBS>::readHofX(ObsEnsemble_ & obsens, const size_t nens,
                                               const size_t iteration,
                                               Observations_ & y_mean_xb) {
  for (size_t jj = 0; jj < nens; ++jj) {
    obsens[jj].read("hofx"+std::to_string(iteration)+"_"+std::to_string(jj+1));
    Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << obsens[jj] << std::endl;
  }
  y_mean_xb.read("hofx_y_mean_xb"+std::to_string(iteration));
}

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
Observations<OBS> LocalEnsembleSolver<MODEL, OBS>::computeHofX(
                                                   const StateEnsemble4D_ & ens_xx,
                                                   size_t iteration,
                                                   bool readFromDisk) {
  util::Timer timer(classname(), "computeHofX");

  const size_t nens = ens_xx.size();
  ObsEnsemble_ obsens(obspaces_, nens);
  Observations_ y_mean_xb(obspaces_);

  // Initialize R_ anew for each iteration
  R_.reset(new ObsErrors_(this->observersconf_, this->obspaces_));

  if (readFromDisk) {
    readHofX(obsens, nens, iteration, y_mean_xb);
    initializeAssimilatedMask();
  } else {
    // compute and save H(x)

    std::vector<util::DateTime> times =
        useLinearObserver() ? ens_xx[0].validTimes() : xbmean_.validTimes();
    util::Duration flength = times[times.size() - 1] - times[0];
    // default_tstep = 2*observation window is passed to PseudoModel as the default
    // pseudomodel time step. It is only used when StateSet has a single state, to enable
    // processing of all observations in the specified window regardless of where in
    // the time window the state is. Observations in
    // ( max(winbgn, xx.time - tstep/2); min(winend, xx.time + tstep/2) ] are
    // processed in H(x).
    util::Duration default_tstep = (obspaces_.windowEnd() - obspaces_.windowStart()) * 2;
    const ModelAux_ moderr(geometry_, eckit::LocalConfiguration());
    const ModelAuxInc_  moderrinc(geometry_, eckit::LocalConfiguration());
    const ObsAux_  obsaux(obspaces_, observersconf_);
    const ObsAuxInc_  obsauxinc(obspaces_, observersconf_);

    // set up postprocessors for the linear model run on ensemble perturbations
    PostProcessor<Increment_> posttl;
    PostProcessorTLAD<MODEL> posttrajtl;

    // save QC filters, obs bias, ob errors for H(mean(Xb) (H(x_i) are saved separately)
    eckit::LocalConfiguration config;
    config.set("save hofx", false);
    config.set("save qc", true);
    config.set("save obs errors", true);
    config.set("save obs bias", true);
    config.set("iteration", std::to_string(iteration));

    // run the model forecast on the background ensemble mean and optionally linearise the observer
    // about this trajectory. save obs errors and qc flags for H(mean(Xb))
    computeHofX4D(config, xbmean_, y_mean_xb, flength, default_tstep, obsaux, moderr,
                  *R_, qcflags_);
    y_mean_xb.save("hofx_y_mean_xb"+std::to_string(iteration));
    // QC flags and Obs errors are set to that of the H(mean(Xb))
    R_->save("ObsError");
    initializeAssimilatedMask();

    // save hofx means that hofx will be written out into ObsSpace;
    // if run computeHofX4D several times with save hofx on,
    // the hofx will be overwritten,
    // unless each time specifying iteration differently in the passed config.
    config.set("save hofx", false);
    config.set("save qc", false);
    config.set("save obs errors", false);
    config.set("save obs bias", false);

    if (useLinearObserver()) {
      // initialize Yb_ and obsloc_
      Yb_ = std::make_unique<DeparturesEnsemble_>(obspaces_, nens);
      obsloc_ = std::make_unique<ObsLocalizations_>(observersconf_, obspaces_);
      // add linearized H(x) to the linear model postprocessor
      linear_hofx_->initializeTL(posttrajtl);
    }
    Departures_ tmpDeps(this->obspaces_);

    // use temporary objects for QC flags and obs errors for ensemble members
    // to avoid overwriting the ones from the H(mean(Xb)) calculation
    std::vector<ObsDataInt_> qcflags;
    for (size_t jobs = 0; jobs < obspaces_.size(); ++jobs) {
      ObsDataInt_ flags(obspaces_[jobs], obspaces_[jobs].obsvariables());
      qcflags.push_back(flags);
    }
    ObsErrors_ Rmat(observersconf_, obspaces_);

    for (size_t jj = 0; jj < nens; ++jj) {
      if (useLinearObserver()) {
        // Setup PseudoLinearModelIncrement4D to run on ensemble perturbation
        Increment4D_ dx(geometry_, ens_xx[jj].variables(), times);
        dx.diff(ens_xx[jj], xbmean_);

        // Approximate H(x_i) (obsens[jj]) around the ensemble mean using linearized model and
        // linearized observer. Firstly, apply the linearized obs operator to this ensemble member's
        // background perturbation from the ensemble mean.
        applyLinearToPerturbations(dx, flength, default_tstep, obsauxinc, moderrinc, posttl,
                                   posttrajtl, tmpDeps);
        // Secondly, add this to the ensemble mean in observation space (calculated with the
        // nonlinear obs operator) giving the approximate H(x_i)
        Yb_->setData(jj, tmpDeps);
        obsens[jj] = y_mean_xb;
        obsens[jj] += Yb_->getData(jj);
      } else {
        // These are recalculated for each ensemble member
        times = ens_xx[jj].validTimes();
        flength = times[times.size()-1] - times[0];
        default_tstep = (obspaces_.windowEnd() - obspaces_.windowStart()) * 2;

        const ModelAux_ moderr(geometry_, eckit::LocalConfiguration());
        const ModelAuxInc_  moderrinc(geometry_, eckit::LocalConfiguration());
        const ObsAux_  obsaux(obspaces_, observersconf_);
        const ObsAuxInc_  obsauxinc(obspaces_, observersconf_);
        computeHofX4D(config, ens_xx[jj], obsens[jj], flength, default_tstep, obsaux, moderr,
                      Rmat, qcflags);
      }
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << obsens[jj] << std::endl;
      obsens[jj].save("hofx"+std::to_string(iteration)+"_"+std::to_string(jj+1));
    }
  }
  // calculate H(x) ensemble mean
  Observations_ yb_mean(obsens.mean());

  // treat the special case of nens=1
  // default option: xbmean_=mean(xb) then yb_mean == y_mean_xb and action below is a tautology
  // if use control member==true: xbmean_ was read from the controll member,
  //                              then using H(xbmean_) is expected by downstream applications
  if (nens == 1) {yb_mean = y_mean_xb;}

  // calculate obs departures
  Observations_ yobs(obspaces_, "ObsValue");
  omb_ = yobs - yb_mean;
  // Need to mask out any missing departures as well as those that have failed QC
  updateAssimilatedMask(omb_);

  // Calculating observation ensemble perturbations when either HofX is read in or a non-linear
  // observer is used. When a linear observer is used, observation ensemble perturbations are
  // implicitly calculated.
  // Also make sure that obs that have missing values in one ensemble member fail for all
  // (this is for the case where different QC procedures are done on different
  // ensemble members).
  Departures_ tmpDeps(this->obspaces_);
  if (readFromDisk || (!useLinearObserver())) {
    // initialize Yb_ and obsloc_
    Yb_ = std::make_unique<DeparturesEnsemble_>(obspaces_, nens);
    obsloc_ = std::make_unique<ObsLocalizations_>(observersconf_, obspaces_);
    for (size_t iens = 0; iens < Yb_->size(); ++iens) {
      tmpDeps = obsens[iens] - yb_mean;
      updateAssimilatedMask(tmpDeps);
      Yb_->setData(iens, tmpDeps);
    }
  }
  // apply assimilated mask to the observation ensemble perturbations
  for (size_t iens = 0; iens < Yb_->size(); ++iens) {
    tmpDeps = Yb_->getData(iens);
    applyAssimilatedMask(tmpDeps);
    Yb_->setData(iens, tmpDeps);
  }
  // apply assimilated mask to the mean departures
  applyAssimilatedMask(omb_);

  // return mean H(x)
  return yb_mean;
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
