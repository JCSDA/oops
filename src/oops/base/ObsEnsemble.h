/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSENSEMBLE_H_
#define OOPS_BASE_OBSENSEMBLE_H_

#include <vector>

#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Ensemble of observations (can hold ensemble of H(x))
template<typename OBS> class ObsEnsemble {
  typedef Observations<OBS>        Observations_;
  typedef ObsSpaces<OBS>           ObsSpaces_;

 public:
  /// Create ensemble of empty Observations size \p nens
  ObsEnsemble(const ObsSpaces_ &, const size_t & nens);

  /// Accessors and size
  size_t size() const {return ensemble_.size();}
  Observations_ & operator[](const size_t ii) {return ensemble_[ii];}
  const Observations_ & operator[](const size_t ii) const { return ensemble_[ii];}

  /// Compute ensemble mean
  Observations_ mean() const;

 private:
  const ObsSpaces_ & obsdb_;             // ObsSpaces used for creating ensemble members
  std::vector<Observations_> ensemble_;  // ensemble members
};

// ====================================================================================

template<typename OBS>
ObsEnsemble<OBS>::ObsEnsemble(const ObsSpaces_ & obsdb, const size_t & nens)
  : obsdb_(obsdb), ensemble_()
{
  ensemble_.reserve(nens);
  for (size_t iens = 0; iens < nens ; ++iens) {
    ensemble_.emplace_back(obsdb_);
  }
  Log::trace() << "ObsEnsemble created" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
Observations<OBS> ObsEnsemble<OBS>::mean() const {
  Observations_ mean_obs(obsdb_);
  mean_obs.zero();
  for (const auto & yy : ensemble_) {
    mean_obs.accumul(yy);
  }
  mean_obs *= 1.0/static_cast<float>(ensemble_.size());
  Log::trace() << "ObsEnsemble::mean done" << std::endl;
  return mean_obs;
}

}  // namespace oops

#endif  // OOPS_BASE_OBSENSEMBLE_H_
