/*
 * (C) Copyright 2018-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_PSEUDOMODEL_H_
#define OOPS_GENERIC_PSEUDOMODEL_H_

#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/util/ConfigFunctions.h"

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/ModelBase.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Configuration handling for ensemble states from template.
class MultipleStateTemplateConfigValidator {
 public:
  explicit MultipleStateTemplateConfigValidator(eckit::LocalConfiguration stateTemplateConf) {
    // template to define a generic state
    config.set("template", stateTemplateConf.getSubConfiguration("template"));
    // Starting state datetime for the first state
    config.set("start datetime", stateTemplateConf.getString("start datetime"));
    // total number of states in the list after excluding certain values
    config.set("number of states", stateTemplateConf.getUnsigned("number of states"));

    // pattern to be replaced by a padded integer for each state
    if (stateTemplateConf.has("pattern")) {
      config.set("pattern", stateTemplateConf.getString("pattern"));
    }

    // starting state padded integer value for the first state
    config.set("start", size_t(1));
    if (stateTemplateConf.has("start")) {
      config.set("start", stateTemplateConf.getUnsigned("start"));
    }
    // list of excluded state pattern values
    config.set("except", std::vector<size_t>{});
    if (stateTemplateConf.has("except")) {
      config.set("except", stateTemplateConf.getUnsignedVector("except"));
    }

    // zero padding
    config.set("zero padding", size_t(0));
    if (stateTemplateConf.has("zero padding")) {
      config.set("zero padding", stateTemplateConf.getUnsigned("zero padding"));
    }
  }

  /// \brief Get State parameters for a given index
  eckit::LocalConfiguration indexToStateConfig(const size_t & stateConfIndex,
      const util::Duration tstep) const {
    // Get correct index
    size_t patternValue = this->config.getUnsigned("start");
    util::DateTime dateTimeValue(this->config.getString("start datetime"));
    std::vector<size_t> except(this->config.getUnsignedVector("except"));
    for (size_t jj = 0; jj <= stateConfIndex; ++jj) {
      // Check for excluded states
      while (std::count(except.begin(), except.end(), patternValue)) {
        patternValue += 1;
      }
      // Update patternValue
      if (jj < stateConfIndex) {
        patternValue += 1;
        dateTimeValue += tstep;
      }
    }

    // Copy and update template configuration with the current pattern value and date
    eckit::LocalConfiguration stateConf(this->config.getSubConfiguration("template"));
    stateConf.set("date", dateTimeValue.toString());
    if (this->config.has("pattern")) {
      std::string pattern(this->config.getString("pattern"));
      util::seekAndReplace(stateConf, pattern, patternValue,
                           this->config.getUnsigned("zero padding"));
    }

    return stateConf;
  }

 private:
  eckit::LocalConfiguration config;
};

// -----------------------------------------------------------------------------
/// \brief Configuration handling for  a list of states.
class MultipleStateConfigValidator {
 public:
  explicit MultipleStateConfigValidator(eckit::LocalConfiguration multipleStateConf) {
    // list of the states in the grouping
    if (multipleStateConf.has("states")) {
      config.set("states", multipleStateConf.getSubConfigurations("states"));
    }
    // template to define states of the list of multiple states
    if (multipleStateConf.has("states from template")) {
      config.set("states from template",
          multipleStateConf.getSubConfiguration("states from template"));
    }
    // run identifier
    config.set("ID", 0);
    if (multipleStateConf.has("ID")) {
      config.set("ID", multipleStateConf.getUnsigned("ID"));
    }

    if (!config.has("states") && !config.has("states from template")) {
      throw eckit::UserError(
          "MultipleStateConfigValidator: both 'states' and 'states from template' are missing",
          Here());
    }
    if (config.has("states") && config.has("states from template")) {
      throw eckit::UserError(
          "MultipleStateConfigValidator: both 'states' and 'states from template' are present",
          Here());
    }
  }

  /// Get number of states in the parameter list
  size_t size() const
  {
    if (this->config.has("states")) {
      return this->config.getSubConfigurations("states").size();
    }
    auto templateConfig = this->config.getSubConfiguration("states from template");
    return templateConfig.getUnsigned("number of states");
  }

  /// Get State parameters for a given index
  eckit::LocalConfiguration getStateConfig(const size_t & stateConfIndex,
      const util::Duration tstep) const {
    // Check number of states
    if (stateConfIndex >= this->size()) {
      throw eckit::UserError(
                "MultipleStateConfigValidator::getStateConfig index argument is too large",
                Here());
    }

    if (this->config.has("states")) {
      // Explicit states
      auto statesConfigs = this->config.getSubConfigurations("states");
      return statesConfigs[stateConfIndex];
    }
    MultipleStateTemplateConfigValidator validator{
      this->config.getSubConfiguration("states from template")
    };
    return validator.indexToStateConfig(stateConfIndex, tstep);
  }

  const eckit::LocalConfiguration & get() const {return this->config;}

 private:
  eckit::LocalConfiguration config;
};

// -----------------------------------------------------------------------------
/// \brief Configuration handling for a pseudomodel (contains either a list of states,
///        or multiple runs).
class PseudoModelConfigValidator {
 public:
  explicit PseudoModelConfigValidator(eckit::LocalConfiguration pseudoModelConf) {
    // list of the states in the grouping
    if (pseudoModelConf.has("states")) {
      config.set("states", pseudoModelConf.getSubConfigurations("states"));
    }
    // template to define states of the list of multiple states
    if (pseudoModelConf.has("states from template")) {
      config.set("states from template",
        pseudoModelConf.getSubConfiguration("states from template"));
    }
    // list of the runs in the Pseudo-model
    if (pseudoModelConf.has("multiple runs")) {
      config.set("multiple runs",
        pseudoModelConf.getSubConfigurations("multiple runs"));
    }

    // Name of the model type
    if (pseudoModelConf.getString("name") != "PseudoModel") {
        throw eckit::UserError(
            "PseudoModelConfigValidator: 'name' is not 'PseudoModel' for this PseudoModel"
            " configuration.",
            Here());
    }
    config.set("name", "PseudoModel");
    // time step between the states in a run
    config.set("tstep", pseudoModelConf.getString("tstep"));

    if (!config.has("states")
        && !config.has("states from template")
        && !config.has("multiple runs")) {
      throw eckit::UserError(
          "PseudoModelConfigValidator: 'states', 'states from template', and 'multiple runs'"
          " parameters are all missing. At least one must be specified.",
          Here());
    }
    if (config.has("states") && config.has("states from template")) {
      throw eckit::UserError(
          "PseudoModelConfigValidator: both states and states from template parameters are"
          " present. Only one or the other should be specified.",
          Here());
    }
  }

  const eckit::LocalConfiguration & get() const {return this->config;}

 private:
  eckit::LocalConfiguration config;
};

// -----------------------------------------------------------------------------

///  Generic implementation of the pseudo model (steps through time by reading states)
template <typename MODEL>
class PseudoModel : public ModelBase<MODEL> {
  typedef Geometry<MODEL>          Geometry_;
  typedef ModelAuxControl<MODEL>   ModelAux_;
  typedef State<MODEL>             State_;

 public:
  static const std::string classname() {return "oops::PseudoModel";}

  PseudoModel(const Geometry_ & resol, const eckit::Configuration & config);

/// initialize forecast
  void initialize(State_ &) const override;
/// one forecast step
  void step(State_ &, const ModelAux_ &) const override;
/// finalize forecast
  void finalize(State_ &) const override;

/// model time step
  const util::Duration & timeResolution() const override {return tstep_;}

 private:
  void print(std::ostream &) const override;
  const util::Duration tstep_;
  std::vector<std::vector<eckit::LocalConfiguration>> runs_;
  mutable std::vector<size_t> currentState_;
  mutable std::vector<util::DateTime> previousTime_;
  mutable size_t currentID_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
PseudoModel<MODEL>::PseudoModel(const Geometry_ & resol, const eckit::Configuration & config)
: tstep_(config.getString("tstep")) {
  PseudoModelConfigValidator validator(eckit::LocalConfiguration{config});
  auto validatedConfig = validator.get();
  if (validatedConfig.has("multiple runs")) {
    // add multiple runs
    auto runConfigs = validatedConfig.getSubConfigurations("multiple runs");
    size_t IDChecker(0);
    for (const auto & runConfig : runConfigs) {
      MultipleStateConfigValidator runValidator(runConfig);
      auto validatedRunConfig = runValidator.get();
      ASSERT(validatedRunConfig.getUnsigned("ID") == IDChecker++);
      std::vector<eckit::LocalConfiguration> runStates;
      for (size_t iState = 0; iState < runValidator.size(); ++iState) {
        runStates.emplace_back(runValidator.getStateConfig(iState, tstep_));
      }
      runs_.emplace_back(runStates);
      currentState_.push_back(0);
      previousTime_.push_back(util::DateTime());
    }
  } else if (validatedConfig.has("states") ||
             validatedConfig.has("states from template")) {
    // add only a single list of states
    MultipleStateConfigValidator runValidator(eckit::LocalConfiguration{config});
    std::vector<eckit::LocalConfiguration> runStates;
    for (size_t iState = 0; iState < runValidator.size(); ++iState) {
      runStates.emplace_back(runValidator.getStateConfig(iState, tstep_));
    }
    runs_.emplace_back(runStates);
    currentState_.push_back(0);
    previousTime_.push_back(util::DateTime());
  }
  Log::trace() << "PseudoModel<MODEL>::PseudoModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::initialize(State_ & xx) const {
  Log::trace() << "PseudoModel<MODEL>::initialize starting" << std::endl;
  currentID_ = xx.ID();

  // If forecast is starting from initial condition, set previousTime_ accordingly
  if (!previousTime_[currentID_].isSet()) previousTime_[currentID_] = xx.validTime();

  // Check whether forecast is being restarted from initial condition before finishing
  if (xx.validTime() < previousTime_[currentID_]) {
    // If it is, reinitialise previousTime_ and currentState_
    previousTime_[currentID_] = xx.validTime();
    currentState_[currentID_] = 0;
  }
  Log::trace() << "PseudoModel<MODEL>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::step(State_ & xx, const ModelAux_ & merr) const {
  Log::trace() << "PseudoModel<MODEL>::step starting" << std::endl;
  // Check that number of steps called for has not exceeded number of State
  // configurations available
  if (currentState_[currentID_] >= runs_[currentID_].size()) {
    std::ostringstream msg;
    msg << classname() << "::step: number of steps called for (" << currentState_[currentID_]
        << ") has exceeded number of State configurations available ("
        << runs_[currentID_].size() << ")" << std::endl;
    throw eckit::UserError(msg.str(), Here());
  }

  xx.updateTime(tstep_);
  xx.read(runs_[currentID_][currentState_[currentID_]++]);  // currentState_ post-incremented

  // Check that time difference between previous State and current State matches tstep_
  // (i.e. that configuration passed to State::read() was for a State at the correct time)
  if (xx.validTime() - previousTime_[currentID_] != tstep_) {
    std::ostringstream msg;
    msg << classname() << "::step: time difference between previous state ("
        << previousTime_[currentID_] << ") and current state (" << xx.validTime()
        << ") does not match tstep_ (" << tstep_ << ")" << std::endl;
    throw eckit::UserError(msg.str(), Here());
  }

  previousTime_[currentID_] = xx.validTime();
  Log::trace() << "PseudoModel<MODEL>::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::finalize(State_ & xx) const {
  Log::trace() << "PseudoModel<MODEL>::finalize starting" << std::endl;
  // If forecast has reached final State configuration, reinitialise variables
  // ready for reforecast
  if (currentState_[currentID_] == runs_[currentID_].size()) {
    currentState_[currentID_] = 0;
    previousTime_[currentID_] = util::DateTime();
  }
  Log::trace() << "PseudoModel<MODEL>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::print(std::ostream & os) const {
  os << "Pseudo model reading states from files with " << tstep_ << " time resolution";
}

}  // namespace oops

#endif  // OOPS_GENERIC_PSEUDOMODEL_H_
