/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_STATEENSEMBLE_H_
#define OOPS_BASE_STATEENSEMBLE_H_

#include <string>
#include <utility>
#include <vector>

#include "oops/base/Accumulator.h"
#include "oops/base/State.h"
#include "oops/base/StateParametersND.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

template<typename MODEL> class Geometry;
class Variables;

// -----------------------------------------------------------------------------
/// Parameters for ensemble members from template.
template <typename MODEL>
class StateMemberTemplateParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(StateMemberTemplateParameters, Parameters)

  typedef StateParametersND<MODEL> Parameters_;
 public:
  RequiredParameter<Parameters_> state{"template", "template to define a generic member", this};
  RequiredParameter<std::string> pattern{"pattern", "pattern to be replaced for members", this};
  RequiredParameter<size_t> nmembers{"nmembers", "number of members", this};
  Parameter<size_t> start{"start", "starting member index", 1, this};
  Parameter<std::vector<size_t>> except{"except", "excluded members indices", {}, this};
  Parameter<size_t> zpad{"zero padding", "zero padding", 0, this};
};

/// Parameters for the ensemble of states.
template <typename MODEL>
class StateEnsembleParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(StateEnsembleParameters, Parameters)

  typedef StateParametersND<MODEL> Parameters_;
  typedef StateMemberTemplateParameters<MODEL> StateMemberTemplateParameters_;
 public:
  OptionalParameter<std::vector<Parameters_>> states{"members",
                   "members of the state ensemble", this};
  OptionalParameter<StateMemberTemplateParameters_> states_template{"members from template",
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

  /// Calculate ensemble mean
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
  // Abort if both "members" and "members from template" are specified
  if (params.states.value() != boost::none && params.states_template.value() != boost::none)
    ABORT("StateEnsemble:contructor: both members and members from template are specified");

  if (params.states.value() != boost::none) {
    // Explicit members

    // Reserve memory to hold ensemble
    const size_t nens = params.states.value()->size();
    states_.reserve(nens);

    // Loop over all ensemble members
    for (size_t jj = 0; jj < nens; ++jj) {
      states_.emplace_back(State_(resol, (*params.states.value())[jj]));
    }
  } else if (params.states_template.value() != boost::none) {
    // Members template

    // Template configuration
    eckit::LocalConfiguration stateConf;
    params.states_template.value()->state.value().serialize(stateConf);

    // Reserve memory to hold ensemble
    const size_t nens = params.states_template.value()->nmembers.value();
    states_.reserve(nens);

    // Loop over all ensemble members
    size_t count = params.states_template.value()->start;
    for (size_t jj = 0; jj < nens; ++jj) {
      // Check for excluded members
      while (std::count(params.states_template.value()->except.value().begin(),
             params.states_template.value()->except.value().end(), count)) {
        count += 1;
      }

      // Copy and update template configuration with pattern
      eckit::LocalConfiguration memberConf(stateConf);

      // Replace pattern recursively in the configuration
      util::seekAndReplace(memberConf, params.states_template.value()->pattern,
        count, params.states_template.value()->zpad);

      // Read state
      states_.emplace_back(State_(resol, memberConf));

      // Update counter
      count += 1;
    }
  } else {
    ABORT("StateEnsemble:contructor: ensemble not specified");
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
