/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSERVERS_H_
#define OOPS_BASE_OBSERVERS_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/GetValuePosts.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observer.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ObsVector.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Computes observation operator (from GeoVaLs), applies bias correction
///        and runs QC filters
template <typename MODEL, typename OBS>
class Observers {
  typedef Geometry<MODEL>               Geometry_;
  typedef GetValuePosts<MODEL, OBS>     GetValuePosts_;
  typedef ObsAuxControls<OBS>           ObsAuxCtrls_;
  typedef ObsErrors<OBS>                ObsErrors_;
  typedef Observations<OBS>             Observations_;
  typedef Observer<MODEL, OBS>          Observer_;
  typedef ObsSpaces<OBS>                ObsSpaces_;
  typedef ObsVector<OBS>                ObsVector_;
  typedef State<MODEL>                  State_;
  typedef PostProcessor<State_>         PostProc_;

 public:
/// \brief Initializes ObsOperators, Locations, and QC data
  Observers(const ObsSpaces_ &, const eckit::Configuration &);

/// \brief Initializes variables, obs bias, obs filters (could be different for
/// different iterations
  void initialize(const Geometry_ &, const ObsAuxCtrls_ &, ObsErrors_ &,
                  PostProc_ &, const int iter = -1);

/// \brief Computes H(x) from the filled in GeoVaLs
  void finalize(Observations_ &);

 private:
  std::vector<std::unique_ptr<Observer_>>  observers_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observers<MODEL, OBS>::Observers(const ObsSpaces_ & obspaces, const eckit::Configuration & config)
  : observers_()
{
  Log::trace() << "Observers<MODEL, OBS>::Observers start" << std::endl;

  std::vector<eckit::LocalConfiguration> obsconfs = config.getSubConfigurations();
  for (size_t jj = 0; jj < obspaces.size(); ++jj) {
    observers_.emplace_back(new Observer_(obspaces[jj], obsconfs[jj]));
  }

  Log::trace() << "Observers<MODEL, OBS>::Observers done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::initialize(const Geometry_ & geom, const ObsAuxCtrls_ & obsaux,
                                       ObsErrors_ & Rmat, PostProc_ & pp, const int iter) {
  Log::trace() << "Observers<MODEL, OBS>::initialize start" << std::endl;

  std::shared_ptr<GetValuePosts_> getvals(new GetValuePosts_());
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    getvals->append(observers_[jj]->initialize(geom, obsaux[jj], Rmat[jj], iter));
  }
  pp.enrollProcessor(getvals);

  Log::trace() << "Observers<MODEL, OBS>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::finalize(Observations_ & yobs) {
  oops::Log::trace() << "Observers<MODEL, OBS>::finalize start" << std::endl;

  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    observers_[jj]->finalize(yobs[jj]);
  }

  oops::Log::trace() << "Observers<MODEL, OBS>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERS_H_
