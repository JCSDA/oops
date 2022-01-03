/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_STATEENSEMBLE_H_
#define OOPS_BASE_STATEENSEMBLE_H_

#include <utility>
#include <vector>

#include "oops/base/Accumulator.h"
#include "oops/base/State.h"
#include "oops/base/StateParametersND.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

template<typename MODEL> class Geometry;
class Variables;

// -----------------------------------------------------------------------------
/// Parameters for the ensemble of states.
template <typename MODEL>
class StateEnsembleParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(StateEnsembleParameters, Parameters)

  typedef StateParametersND<MODEL> Parameters_;
 public:
  RequiredParameter<std::vector<Parameters_>> states{"members",
                   "members of the state ensemble", this};
};

/// \brief Ensemble of states
template<typename MODEL> class StateEnsemble {
  typedef Geometry<MODEL>      Geometry_;
  typedef State<MODEL>         State_;
  typedef StateEnsembleParameters<MODEL> StateEnsembleParameters_;

 public:
  /// Create ensemble of states
  StateEnsemble(const Geometry_ &, const StateEnsembleParameters_ &);

  /// calculate ensemble mean
  State_ mean() const;

  /// Accessors
  size_t size() const { return states_.size(); }
  State_ & operator[](const int ii) { return states_[ii]; }
  const State_ & operator[](const int ii) const { return states_[ii]; }

  /// Information
  const Variables & variables() const {return states_[0].variables();}

 private:
  std::vector<State_> states_;
};

// ====================================================================================

template<typename MODEL>
StateEnsemble<MODEL>::StateEnsemble(const Geometry_ & resol,
                                    const StateEnsembleParameters_ & params)
  : states_() {
  const size_t nens = params.states.value().size();
  states_.reserve(nens);
  // Loop over all ensemble members
  for (size_t jj = 0; jj < nens; ++jj) {
    states_.emplace_back(State_(resol, params.states.value()[jj]));
  }
  Log::trace() << "StateEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> StateEnsemble<MODEL>::mean() const {
  // Compute ensemble mean
  Accumulator<MODEL, State_, State_> ensmean(states_[0]);

  const double rr = 1.0/static_cast<double>(states_.size());
  for (size_t iens = 0; iens < states_.size(); ++iens) {
    ensmean.accumul(rr, states_[iens]);
  }

  Log::trace() << "StateEnsemble::mean done" << std::endl;
  return std::move(ensmean);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_STATEENSEMBLE_H_
