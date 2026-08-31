/*
 * (C) Copyright 2020 UCAR.
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSERVERS_H_
#define OOPS_BASE_OBSERVERS_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Geometry.h"
#include "oops/base/GetValuePerts.h"
#include "oops/base/GetValuePosts.h"
#include "oops/base/GetValueTLADs.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observer.h"
#include "oops/base/ObsOperatorBase.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/ObsVector.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Computes observation operator (from GeoVaLs), applies bias correction
///        and runs QC filters
template <typename MODEL, typename OBS>
class Observers {
  typedef Geometry<MODEL>               Geometry_;
  typedef GetValues<MODEL, OBS>         GetValues_;
  typedef GetValuePerts<MODEL, OBS>     GetValuePerts_;
  typedef GetValuePosts<MODEL, OBS>     GetValuePosts_;
  typedef GetValueTLADs<MODEL, OBS>     GetValueTLADs_;
  typedef ObsAuxControls<OBS>           ObsAuxCtrls_;
  typedef ObsDataVector<OBS, int>       ObsDataInt_;
  typedef ObsErrors<OBS>                ObsErrors_;
  typedef Observations<OBS>             Observations_;
  typedef Observer<MODEL, OBS>          Observer_;
  typedef ObsOperatorBase<OBS>          ObsOperatorBase_;
  typedef ObsSpaces<OBS>                ObsSpaces_;
  typedef ObsVector<OBS>                ObsVector_;
  typedef State<MODEL>                  State_;
  typedef PostProcessor<State_>         PostProc_;
  template <typename DATA> using ObsData_ = ObsDataVector<OBS, DATA>;
  template <typename DATA> using ObsDataVec_ = std::vector<ObsData_<DATA>>;

 public:
  Observers(const ObsSpaces_ &, const eckit::Configuration &,
            std::vector<std::unique_ptr<ObsOperatorBase_>> obsOpBases = {});

/// \brief Initializes variables, obs bias, obs filters (could be different for
/// different iterations
  void initialize(const Geometry_ &, ObsAuxCtrls_ &, ObsErrors_ &,
                  PostProc_ &, const eckit::Configuration & = eckit::LocalConfiguration());

/// \brief Computes H(x) from the filled in GeoVaLs
  void finalize(Observations_ &, std::vector<ObsDataInt_> &);

  void resetObsPert(const Geometry_ &, std::vector<std::unique_ptr<ObsOperatorBase_>>,
                    const std::shared_ptr<GetValueTLADs_> &, const Variables &);

  void updateObservers(const eckit::Configuration &);

 private:
  std::vector<std::unique_ptr<Observer_>>  observers_;
  eckit::LocalConfiguration getValuesConf_;
  std::vector<eckit::LocalConfiguration> obsconfs_;
  std::shared_ptr<GetValuePosts_> posts_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observers<MODEL, OBS>::Observers(const ObsSpaces_ & obspaces, const eckit::Configuration & config,
                                 std::vector<std::unique_ptr<ObsOperatorBase_>> obsOpBases)
  : observers_(), getValuesConf_(config.getSubConfiguration("get values")),
    obsconfs_(config.getSubConfigurations("observers")),
    posts_(new GetValuePosts_(getValuesConf_))
{
  Log::trace() << "Observers<MODEL, OBS>::Observers start" << std::endl;

  if (obsOpBases.size() != obspaces.size()) obsOpBases.resize(obspaces.size());
  ASSERT(obspaces.size() == obsconfs_.size());
  for (size_t jj = 0; jj < obspaces.size(); ++jj) {
    observers_.emplace_back(new Observer_(obspaces[jj], obsconfs_[jj], std::move(obsOpBases[jj])));
  }

  Log::trace() << "Observers<MODEL, OBS>::Observers done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::initialize(const Geometry_ & geom, ObsAuxCtrls_ & obsaux,
                                       ObsErrors_ & Rmat, PostProc_ & pp,
                                       const eckit::Configuration & conf) {
  Log::trace() << "Observers<MODEL, OBS>::initialize start" << std::endl;

  posts_->clear();
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    auto getval = observers_[jj]->initialize(geom, obsaux[jj], Rmat[jj], conf);
    if (getval) {
      // getval will be nullptr if observer is initialized with a GeoVaLs file
      posts_->append(getval);
    }
  }
  pp.enrollProcessor(posts_);

  Log::trace() << "Observers<MODEL, OBS>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::finalize(Observations_ & yobs, std::vector<ObsDataInt_> & qcflags) {
  oops::Log::trace() << "Observers<MODEL, OBS>::finalize start" << std::endl;

  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    observers_[jj]->finalize(yobs[jj], qcflags[jj]);
  }

  oops::Log::trace() << "Observers<MODEL, OBS>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::resetObsPert(const Geometry_ & geom,
                                         std::vector<std::unique_ptr<ObsOperatorBase_>> obsOpBases,
                                         const std::shared_ptr<GetValueTLADs_> & getValTLs,
                                         const Variables & vars) {
  oops::Log::trace() << "Observers<MODEL, OBS>::resetObsOp start" << std::endl;

  posts_.reset(new GetValuePerts_(getValuesConf_, getValTLs, vars));
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    observers_[jj]->resetObsPert(geom, std::move(obsOpBases[jj]), (*getValTLs)[jj]);
  }

  oops::Log::trace() << "Observers<MODEL, OBS>::resetObsOp done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::updateObservers(const eckit::Configuration & cdaConfig) {
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    observers_[jj]->updateObserver(cdaConfig);
  }
  Log::trace() << "Observers obs error appended" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERS_H_
