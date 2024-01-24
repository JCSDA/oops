/*
 * (C) Copyright 2019-2020 UCAR.
 * (C) Crown copyright 2024, Met Office
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
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/FieldSetOperations.h"
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

 public:
  RequiredParameter<eckit::LocalConfiguration> state{"template",
                                                     "template to define a generic member", this};
  RequiredParameter<std::string> pattern{"pattern", "pattern to be replaced for members", this};
  RequiredParameter<size_t> nmembers{"nmembers", "number of members", this};
  Parameter<size_t> start{"start", "starting member index", 1, this};
  Parameter<std::vector<size_t>> except{"except", "excluded members indices", {}, this};
  Parameter<size_t> zpad{"zero padding", "zero padding", 0, this};
};

// -----------------------------------------------------------------------------

/// Parameters for the ensemble of states.
template <typename MODEL>
class StateEnsembleParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(StateEnsembleParameters, Parameters)

  typedef StateMemberTemplateParameters<MODEL> StateMemberTemplateParameters_;
 public:
  OptionalParameter<std::vector<eckit::LocalConfiguration>> states{"members",
                   "members of the state ensemble", this};
  OptionalParameter<StateMemberTemplateParameters_> states_template{"members from template",
                   "template to define members of the state ensemble", this};

  /// Overridden to detect missing conditionally required parameters
  using Parameters::deserialize;
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  /// Get ensemble size
  size_t size() const;

  /// Get Increment parameters for a given ensemble index
  eckit::LocalConfiguration getStateConfig(const size_t &, const size_t &) const;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
void StateEnsembleParameters<MODEL>::deserialize(util::CompositePath &path,
                                                 const eckit::Configuration &config)
{
  Parameters::deserialize(path, config);

  if (states.value() == boost::none && states_template.value() == boost::none) {
    throw eckit::UserError(
        path.path() +
        ": both members and members from template are missing",
        Here());
  }
  if (states.value() != boost::none && states_template.value() != boost::none) {
    throw eckit::UserError(
        path.path() +
        ": both members and members from template are present",
        Here());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
size_t StateEnsembleParameters<MODEL>::size() const
{
  if (states.value() != boost::none) {
    return states.value()->size();
  } else {
    return states_template.value()->nmembers.value();
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
eckit::LocalConfiguration StateEnsembleParameters<MODEL>::getStateConfig(const size_t & ie,
                                                                         const size_t & rank) const
{
  // Check ensemble size
  if (ie >= this->size()) {
    ABORT("StateEnsembleParameters: getStateConfig member index is too large");
  }

  if (states.value() != boost::none) {
    // Explicit members
    return (*states.value())[ie];
  } else {
    // Members template

    // Template configuration
    eckit::LocalConfiguration stateConf(states_template.value()->state.value());

    // Get correct index
    size_t count = states_template.value()->start;
    for (size_t jj = 0; jj <= ie; ++jj) {
      // Check for excluded members
      while (std::count(states_template.value()->except.value().begin(),
             states_template.value()->except.value().end(), count)) {
        count += 1;
      }

      // Update counter
      if (jj < ie) count += 1;
    }

    // Copy and update template configuration with pattern
    eckit::LocalConfiguration memberConf(stateConf);

    // Replace pattern recursively in the configuration
    util::seekAndReplace(memberConf, states_template.value()->pattern,
      count, states_template.value()->zpad);

    if (memberConf.has("states")) {
      std::vector<eckit::LocalConfiguration> confs = memberConf.getSubConfigurations("states");
      ASSERT(rank < confs.size());
      memberConf = confs[rank];
    } else {
      ASSERT(rank == 0);
    }

    return memberConf;
  }
}

// -----------------------------------------------------------------------------

/// \brief Ensemble of states
template<typename MODEL> class StateEnsemble {
  typedef Geometry<MODEL>      Geometry_;
  typedef Increment<MODEL>     Increment_;
  typedef State<MODEL>         State_;
  typedef StateEnsembleParameters<MODEL> StateEnsembleParameters_;

 public:
  /// Create ensemble of states
  StateEnsemble(const Geometry_ &, const StateEnsembleParameters_ &);
  /// Create n copies of a state in an ensemble
  StateEnsemble(const State_ &, const size_t &);

  /// Calculate ensemble mean
  State_ mean() const;

  /// Calculate ensemble spread
  Increment_ variance() const;
  Increment_ stddev() const;

  /// Accessors
  size_t size() const { return states_.size(); }
  State_ & operator[](const int ii) { return states_[ii]; }
  const State_ & operator[](const int ii) const { return states_[ii]; }

  /// Information
  const Variables & variables() const {return states_[0].variables();}

 private:
  std::vector<State_> states_;
  const Geometry_ & geom_;
};

// ====================================================================================

template<typename MODEL>
StateEnsemble<MODEL>::StateEnsemble(const Geometry_ & resol,
                                    const StateEnsembleParameters_ & params)
  : states_(), geom_(resol) {
  // Reserve memory to hold ensemble
  const size_t nens = params.size();
  states_.reserve(nens);

  const size_t myrank = resol.timeComm().rank();

  // Loop over all ensemble members
  for (size_t jj = 0; jj < nens; ++jj) {
    states_.emplace_back(resol, params.getStateConfig(jj, myrank));
  }
  Log::trace() << "StateEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateEnsemble<MODEL>::StateEnsemble(const State_ & copyState_, const size_t & ensSize_)
  : states_(), geom_(copyState_.geometry()) {
    states_.reserve(ensSize_);
    for (size_t jj = 0; jj < ensSize_; ++jj) {
    states_.emplace_back(copyState_);
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

template<typename MODEL>
Increment<MODEL> StateEnsemble<MODEL>::variance() const {
  ASSERT(states_.size() > 1);
  // Ensemble mean
  State<MODEL> ensmean = this->mean();

  // Compute ensemble variance
  Increment_ ensVar(geom_, this->variables(), ensmean.validTime());
  ensVar.zero();

  const double rr = 1.0/(static_cast<double>(states_.size()) - 1.0);
  Increment_ pert(ensVar);
  for (size_t iens = 0; iens < states_.size(); ++iens) {
    pert.zero();
    pert.diff(states_[iens], ensmean);
    pert.schur_product_with(pert);
    ensVar.axpy(rr, pert);
  }

  Log::trace() << "StateEnsemble:: variance done" << std::endl;
  return ensVar;
}
// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> StateEnsemble<MODEL>::stddev() const {
  ASSERT(states_.size() > 1);
  // Ensemble variance
  Increment<MODEL> ensStdDev = this->variance();

  // Compute ensemble standard deviation
  ensStdDev.fieldSet().sqrt();
  ensStdDev.synchronizeFields();

  Log::trace() << "StateEnsemble:: standard deviation done" << std::endl;
  return ensStdDev;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_STATEENSEMBLE_H_
