/*
 * (C) Copyright 2018-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_PSEUDOMODEL_H_
#define OOPS_GENERIC_PSEUDOMODEL_H_

#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/ModelBase.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// Configuration options taken by the pseudo model.
template <typename MODEL>
class PseudoModelParameters : public ModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(PseudoModelParameters, ModelParametersBase)

 public:
  typedef typename State<MODEL>::Parameters_ StateParameters_;

  RequiredParameter<util::Duration> tstep{"tstep", "Time step", this};
  RequiredParameter<std::vector<StateParameters_>> states{
    "states",
    "List of configuration options used to initialize the model state at each time step",
    this};
};

// -----------------------------------------------------------------------------

///  Generic implementation of the pseudo model (steps through time by reading states)
template <typename MODEL>
class PseudoModel : public ModelBase<MODEL> {
  typedef Geometry<MODEL>          Geometry_;
  typedef ModelAuxControl<MODEL>   ModelAux_;
  typedef State<MODEL>             State_;
  typedef typename State<MODEL>::Parameters_ StateParameters_;

 public:
  typedef PseudoModelParameters<MODEL> Parameters_;

  static const std::string classname() {return "oops::PseudoModel";}

  PseudoModel(const Geometry_ &, const Parameters_ &);

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
  std::vector<StateParameters_> states_;
  mutable size_t currentstate_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
PseudoModel<MODEL>::PseudoModel(const Geometry_ & resol, const Parameters_ & parameters)
  : tstep_(parameters.tstep),
    states_(parameters.states),
    currentstate_(0) {
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
  xx.read(parametersOrConfiguration<HasParameters_<State<MODEL>>::value>(states_[currentstate_++]));
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
