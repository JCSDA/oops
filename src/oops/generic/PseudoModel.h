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
  std::vector<eckit::LocalConfiguration> states_;
  mutable size_t currentstate_;
  mutable util::DateTime previousTime_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
PseudoModel<MODEL>::PseudoModel(const Geometry_ & resol, const eckit::Configuration & config)
  : tstep_(config.getString("tstep")),
    states_(config.getSubConfigurations("states")), currentstate_(0), previousTime_() {
  Log::trace() << "PseudoModel<MODEL>::PseudoModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::initialize(State_ & xx) const {
  currentstate_ = 0;
  Log::trace() << "PseudoModel<MODEL>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::step(State_ & xx, const ModelAux_ & merr) const {
  Log::trace() << "PseudoModel<MODEL>:step Starting " << std::endl;
  xx.updateTime(tstep_);
  if (currentstate_ >= states_.size()) {
    std::ostringstream msg;
    msg << classname() << "::step mismatch between timestep and number of states: "
        << currentstate_ << "/" << states_.size() << " tstep_ = " << tstep_ << std::endl;
    throw eckit::UserError(msg.str(), Here());
  }
  xx.read(states_[currentstate_++]);

  // Verify different states are separated by the timestep
  if (currentstate_ != 1 && previousTime_ != xx.validTime()) {
    const util::Duration timeStep = timeResolution();
    if (xx.validTime() - previousTime_ != timeStep) {
      std::ostringstream msg;
      msg << classname() << "::step time step does not match time increment between states"
          << " at " << previousTime_ << " and " << xx.validTime()
          << " with time step " << timeStep << std::endl;
      throw eckit::UserError(msg.str(), Here());
    }
  }
  previousTime_ = xx.validTime();

  Log::trace() << "PseudoModel<MODEL>::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::finalize(State_ & xx) const {
  Log::trace() << "PseudoModel<MODEL>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModel<MODEL>::print(std::ostream & os) const {
  os << "Pseudo model reading states from files with " << tstep_ << " time resolution";
}

}  // namespace oops

#endif  // OOPS_GENERIC_PSEUDOMODEL_H_
