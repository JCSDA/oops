/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSENSEMBLE_H_
#define OOPS_BASE_OBSENSEMBLE_H_

#include <memory>
#include <vector>

#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Ensemble

template<typename MODEL> class ObsEnsemble {
  typedef Observations<MODEL>        Observations_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;

 public:
/// Constructor
  ObsEnsemble(const ObsSpaces_ &, const unsigned int &);

/// Destructor
  virtual ~ObsEnsemble() {}

  /// Accessors
  unsigned int size() const {
    return ensemble_.size();
  }
  Observations_ & operator[](const int ii) {
    return *ensemble_[ii];
  }
  const Observations_ & operator[](const int ii) const {
    return *ensemble_[ii];
  }

  Observations_ mean() const;

 private:
  const ObsSpaces_ & obsdb_;
  std::vector<std::shared_ptr<Observations_>> ensemble_;
};

// ====================================================================================

template<typename MODEL>
ObsEnsemble<MODEL>::ObsEnsemble(const ObsSpaces_ & obsdb, const unsigned int & ns)
  : obsdb_(obsdb), ensemble_()
{
  for (unsigned i = 0; i < ns ; ++i) {
    std::shared_ptr<Observations_> y(new Observations_(obsdb_));
    ensemble_.push_back(y);
  }
  Log::trace() << "ObsEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Observations<MODEL> ObsEnsemble<MODEL>::mean() const {
  Observations_ mean_obs(obsdb_);
  mean_obs.zero();
  for (unsigned ii = 0; ii < ensemble_.size(); ++ii) {
    mean_obs.accumul(*ensemble_[ii]);
  }
  mean_obs *= 1.0/static_cast<float>(ensemble_.size());
  return mean_obs;
}

}  // namespace oops

#endif  // OOPS_BASE_OBSENSEMBLE_H_
