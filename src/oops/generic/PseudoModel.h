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

///  Generic implementation of the pseudo model (steps through time by reading states)
template <typename MODEL>
class PseudoModel : public ModelBase<MODEL> {
  typedef Geometry<MODEL>          Geometry_;
  typedef ModelAuxControl<MODEL>   ModelAux_;
  typedef State<MODEL>             State_;

 public:
  static const std::string classname() {return "oops::PseudoModel";}

  PseudoModel(const Geometry_ &, const eckit::Configuration &);

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
  if (config.has("multiple runs")) {
    std::vector<eckit::LocalConfiguration> configs(config.getSubConfigurations("multiple runs"));
    size_t IDChecker(0);
    for (const auto& run : configs) {
      ASSERT(run.getUnsigned("ID") == IDChecker++);
      runs_.push_back(run.getSubConfigurations("states"));
      currentState_.push_back(0);
      previousTime_.push_back(util::DateTime());
    }
  } else if (config.has("states")) {
    runs_.push_back(config.getSubConfigurations("states"));
    currentState_.push_back(0);
    previousTime_.push_back(util::DateTime());
  } else {
    ABORT("PseudoModel<MODEL>::PseudoModel: unsupported configuration");
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
  // Check that number of steps called for has not exceeded number of State configurations available
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
  // If forecast has reached final State configuration, reinitialise variables ready for reforecast
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
