/*
 * (C) Copyright 2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_STATEENSEMBLE_H_
#define OOPS_BASE_STATEENSEMBLE_H_

#include <memory>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Accumulator.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Ensemble
template<typename MODEL> class StateEnsemble {
  typedef Geometry<MODEL>      Geometry_;
  typedef State<MODEL>         State_;

 public:
  /// Constructor
  StateEnsemble(const Geometry_ & resol, const Variables & vars,
                const eckit::Configuration &);

  /// Destructor
  virtual ~StateEnsemble() {}

  /// ensemble mean
  State_ mean() const;

  /// Accessors
  unsigned int size() const {
    return states_.size();
  }
  State_ & operator[](const int ii) {
    return *states_[ii];
  }
  const State_ & operator[](const int ii) const {
    return *states_[ii];
  }

 private:
  std::vector<std::shared_ptr<State_>> states_;
};

// ====================================================================================

template<typename MODEL>
StateEnsemble<MODEL>::StateEnsemble(const Geometry_ & resol,
                                    const Variables & vars,
                                    const eckit::Configuration & config)
  : states_() {
  std::vector<eckit::LocalConfiguration> memberConfig;
  config.get("members", memberConfig);
  // Loop over all ensemble members
  for (size_t jj = 0; jj < memberConfig.size(); ++jj) {
    std::shared_ptr<State_> xx(new State_(resol, vars, memberConfig[jj]));
    states_.push_back(xx);
  }
  Log::trace() << "StateEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> StateEnsemble<MODEL>::mean() const {
  // Compute ensemble mean
  Accumulator<MODEL, State_, State_> ensmean(*states_[0]);

  const double rr = 1.0/static_cast<double>(states_.size());
  for (size_t iens = 0; iens < states_.size(); ++iens) {
    ensmean.accumul(rr, *states_[iens]);
  }

  Log::trace() << "StateEnsemble::mean done" << std::endl;
  return ensmean;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_STATEENSEMBLE_H_
