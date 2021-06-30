/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_STATEENSEMBLE4D_H_
#define OOPS_BASE_STATEENSEMBLE4D_H_

#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Accumulator.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Ensemble of 4D states
template<typename MODEL> class StateEnsemble4D {
  typedef Geometry<MODEL>      Geometry_;
  typedef State4D<MODEL>       State4D_;

 public:
  /// Create ensemble of 4D states
  StateEnsemble4D(const Geometry_ &, const eckit::Configuration &);

  /// calculate ensemble mean
  State4D_ mean() const;

  /// Accessors
  unsigned int size() const { return states_.size(); }
  State4D_ & operator[](const int ii) { return states_[ii]; }
  const State4D_ & operator[](const int ii) const { return states_[ii]; }

  /// Information
  const Variables & variables() const {return states_[0].variables();}

 private:
  std::vector<State4D_> states_;
};

// ====================================================================================

template<typename MODEL>
StateEnsemble4D<MODEL>::StateEnsemble4D(const Geometry_ & resol,
                                        const eckit::Configuration & config)
  : states_() {
  std::vector<eckit::LocalConfiguration> memberConfig;
  config.get("members", memberConfig);
  states_.reserve(memberConfig.size());
  // Loop over all ensemble members
  for (size_t jj = 0; jj < memberConfig.size(); ++jj) {
    states_.emplace_back(State4D_(resol, memberConfig[jj]));
  }
  Log::trace() << "StateEnsemble4D:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State4D<MODEL> StateEnsemble4D<MODEL>::mean() const {
  // Compute ensemble mean
  Accumulator<MODEL, State4D_, State4D_> ensmean(states_[0]);

  const double rr = 1.0/static_cast<double>(states_.size());
  for (size_t iens = 0; iens < states_.size(); ++iens) {
    ensmean.accumul(rr, states_[iens]);
  }

  Log::trace() << "StateEnsemble4D::mean done" << std::endl;
  return std::move(ensmean);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_STATEENSEMBLE4D_H_
