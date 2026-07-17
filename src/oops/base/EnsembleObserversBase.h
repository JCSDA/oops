/*
 * (C) Copyright 2026- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/base/Departures.h"
#include "oops/base/Geometry.h"
#include "oops/base/Model.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/base/Variables.h"
#include "oops/generic/PseudoModelState4D.h"
#include "oops/generic/VerticalLocEV.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// \brief Common infrastructure shared by the ensemble observers used by the local ensemble DA
/// applications (LETKF/GETKF): running the nonlinear observation operator on a StateSet_ (the
/// building block for computing H(x) of the ensemble mean, and, for nonlinear observers, of each
/// ensemble member), and masking/tracking which observations are assimilated.
///
/// Two kinds of ensemble observers derive from this class - \ref LinearEnsembleObservers and
/// \ref NonlinearEnsembleObservers - depending on whether the config option
/// "local ensemble DA.use linear observer" is set. Both of those handle the optional
/// modulated-ensemble (GETKF) case internally, controlled by \ref modulated_, which is set when
/// this class is constructed with a non-null \p xbmean / \p incvars (used by GETKF-family
/// solvers to compute H(x) for the modulated ensemble perturbations used in the vertical
/// localization eigenvector expansion).
template <typename MODEL, typename OBS>
class EnsembleObserversBase {
 protected:
  typedef Departures<OBS>             Departures_;
  typedef Geometry<MODEL>             Geometry_;
  typedef Model<MODEL>                Model_;
  typedef ModelAuxControl<MODEL>      ModelAux_;
  typedef ObsAuxControls<OBS>         ObsAux_;
  typedef ObsDataVector<OBS, int>     ObsDataInt_;
  typedef ObsError<OBS>               ObsError_;
  typedef ObsErrors<OBS>              ObsErrors_;
  typedef Observations<OBS>           Observations_;
  typedef Observers<MODEL, OBS>       Observers_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef PseudoModelState4D<MODEL>   PseudoModel_;
  typedef State<MODEL>                State_;
  typedef StateSet<MODEL>             StateSet_;
  typedef VerticalLocEV<MODEL>        VerticalLocEV_;

 public:
  static const std::string classname() {return "oops::EnsembleObserversBase";}

  /// initialize observers with \p obspaces, \p geometry, full \p config and \p nens ensemble
  /// size. If \p xbmean and \p incvars are provided (non-null), the modulated-ensemble (GETKF)
  /// infrastructure (vertical localization eigenvectors) is set up as well.
  /// If \p priorQcFlags is provided (non-null), observations excluded there (e.g. by an earlier
  /// background ensemble observer) are also excluded from the final ensemble-mean departures
  /// computed by this object (e.g. a later posterior/analysis ensemble observer) - see
  /// \ref maskWithPriorQcFlags. Each ensemble observer otherwise independently (re-)computes its
  /// own QC decisions from scratch, which without this would let an observation excluded from
  /// e.g. the background diagnostics reappear in the posterior ones (or vice versa).
  EnsembleObserversBase(ObsSpaces_ & obspaces, const Geometry_ & geometry,
                        const eckit::Configuration & config, size_t nens,
                        const StateSet_ * xbmean = nullptr, const Variables * incvars = nullptr,
                        const std::vector<ObsDataInt_> * priorQcFlags = nullptr);
  virtual ~EnsembleObserversBase() = default;

  /// quality control flags resulting from this object's computeHofX, one per obs space - meant
  /// to be passed as \p priorQcFlags to a later ensemble observer (see constructor)
  const std::vector<ObsDataInt_> & qcFlags() const {return qcflags_;}

 protected:
  const Geometry_  & geometry_;   ///< Geometry associated with the updated states
  const ObsSpaces_ & obspaces_;   ///< ObsSpaces used in the update
  const eckit::LocalConfiguration obsconf_;  ///< configuration for observations
  const eckit::LocalConfiguration observersconf_;  ///< configuration for observations.observers
  const size_t nens_;              ///< ensemble size

  std::unique_ptr<ObsErrors_>  R_;         ///< observation errors, set in computeHofX method
  std::unique_ptr<Departures_> invVarR_;   ///< inverse observation error variance for assimilated
                                           ///< observations; set in initializeAssimilatedMask
  std::vector<ObsDataInt_> qcflags_;  ///< quality control flags, set by computeHofX

  /// QC flags inherited from another ensemble observer's computeHofX (e.g. background, when this
  /// is the posterior observer) - see constructor and \ref maskWithPriorQcFlags. Empty if none.
  std::vector<ObsDataInt_> priorQcFlags_;

  /// whether the modulated ensemble (GETKF) infrastructure below is active
  const bool modulated_;
  Variables incvars_;                          ///< increment variables (modulated_ only)
  size_t neig_;                                ///< number of vertical localization eigenvectors
  std::unique_ptr<VerticalLocEV_> vertloc_;     ///< vertical localization eigenvectors

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
  /// Mask \p dep (the final ensemble-mean departures) with this object's own QC flags (from
  /// its H(mean(Xb)) computation, \ref qcflags_) and, if this object was constructed with prior
  /// QC flags (see constructor), also with those - so an observation excluded either here or by
  /// an earlier (e.g. background) ensemble observer is consistently excluded from both.
  void maskWithPriorQcFlags(Departures_ & dep) const {
    dep.mask(qcflags_);
    if (!priorQcFlags_.empty()) dep.mask(priorQcFlags_);
  }
  /// Overload for the ensemble-mean H(x) itself (\ref Observations_), printed alongside \p dep
  /// above - kept consistent with it for the same reason.
  void maskWithPriorQcFlags(Observations_ & yy) const {
    yy.mask(qcflags_);
    if (!priorQcFlags_.empty()) yy.mask(priorQcFlags_);
  }

  /// Apply the non-linear observation operator to the state \p xx and store the result in \p yy.
  /// If \p enrollTrajectory is provided, it is called with the internal nonlinear-model
  /// postprocessor before the forecast is run, allowing a caller to enroll additional
  /// postprocessors (e.g. a TrajectorySaver, used by \ref LinearEnsembleObservers to linearize
  /// the observation operator about this trajectory).
  void computeHofX4D(const eckit::Configuration & config, const StateSet_ & xx, Observations_ & yy,
                     const util::Duration & flength, const util::Duration & default_tstep,
                     const ObsAux_ & obsaux, const ModelAux_ & moderr,
                     ObsErrors_ & R, std::vector<ObsDataInt_> & qcflags,
                     const std::function<void(PostProcessor<State_> &)> & enrollTrajectory =
                         std::function<void(PostProcessor<State_> &)>());
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
EnsembleObserversBase<MODEL, OBS>::EnsembleObserversBase(ObsSpaces_ & obspaces,
                                                          const Geometry_ & geometry,
                                                          const eckit::Configuration & config,
                                                          size_t nens,
                                                          const StateSet_ * xbmean,
                                                          const Variables * incvars,
                                                          const std::vector<ObsDataInt_> *
                                                              priorQcFlags)
  : geometry_(geometry),
    obspaces_(obspaces),
    obsconf_(config, "observations"),
    observersconf_(obsconf_, "observers"),
    nens_(nens),
    modulated_(xbmean != nullptr),
    neig_(0) {
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    ObsDataInt_ qcflags(obspaces_[jj], obspaces_[jj].obsvariables());
    qcflags_.push_back(qcflags);
  }
  if (priorQcFlags != nullptr) {
    ASSERT(priorQcFlags->size() == obspaces_.size());
    priorQcFlags_ = *priorQcFlags;
  }
  if (modulated_) {
    ASSERT(incvars != nullptr);
    incvars_ = *incvars;
    vertloc_ = std::make_unique<VerticalLocEV_>(
        config.getSubConfiguration("local ensemble DA.vertical localization"), (*xbmean)[0],
        incvars_);
    neig_ = vertloc_->neig();
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void EnsembleObserversBase<MODEL, OBS>::computeHofX4D(
    const eckit::Configuration & config, const StateSet_ & xx, Observations_ & yy,
    const util::Duration & flength, const util::Duration & default_tstep,
    const ObsAux_ & obsaux, const ModelAux_ & moderr, ObsErrors_ & R,
    std::vector<ObsDataInt_> & qcflags,
    const std::function<void(PostProcessor<State_> &)> & enrollTrajectory) {
  // Setup pseudo model to run on the passed-in state
  State_ init_xx = xx[0];
  std::unique_ptr<PseudoModel_> pseudomodel(new PseudoModel_(xx, default_tstep));
  const Model_ model(std::move(pseudomodel));

  // setup nonlinear postprocessor nonlinear observers
  PostProcessor<State_> post;
  Observers_ hofx(obspaces_, obsconf_);

  // initialize nonlinear model postprocessor
  hofx.initialize(geometry_, obsaux, R, post, config);

  if (enrollTrajectory) {
    enrollTrajectory(post);
  }

  // run nonlinear model
  model.forecast(init_xx, moderr, flength, post);

  // compute nonlinear H(x)
  hofx.finalize(yy, qcflags);
}

// -----------------------------------------------------------------------------

}  // namespace oops
