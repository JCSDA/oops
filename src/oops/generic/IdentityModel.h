/*
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_IDENTITYMODEL_H_
#define OOPS_GENERIC_IDENTITYMODEL_H_

#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/ModelBase.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

class IdentityModelParameters : public ModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(IdentityModelParameters, ModelParametersBase)

 public:
  oops::RequiredParameter<util::Duration> tstep{"tstep", this};
};

/// Generic implementation of identity model
template <typename MODEL>
class IdentityModel : public ModelBase<MODEL> {
  typedef Geometry<MODEL>          Geometry_;
  typedef ModelAuxControl<MODEL>   ModelAux_;
  typedef State<MODEL>             State_;

 public:
  typedef IdentityModelParameters           Parameters_;

  static const std::string classname() {return "oops::IdentityModel";}

  IdentityModel(const Geometry_ &, const IdentityModelParameters &);

/// initialize forecast
  void initialize(State_ &) const override;
/// one forecast step
  void step(State_ &, const ModelAux_ &) const override;
/// finalize forecast
  void finalize(State_ &) const override;

/// model time step
  const util::Duration & timeResolution() const override {return params_.tstep;}

 private:
  void print(std::ostream &) const override {}
  const IdentityModelParameters params_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
IdentityModel<MODEL>::IdentityModel(const Geometry_ & resol, const IdentityModelParameters & params)
  : params_(params) {
  Log::trace() << "IdentityModel<MODEL>::IdentityModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityModel<MODEL>::initialize(State_ & xx) const {
  Log::trace() << "IdentityModel<MODEL>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityModel<MODEL>::step(State_ & xx, const ModelAux_ & merr) const {
  Log::trace() << "IdentityModel<MODEL>:step Starting " << std::endl;
  xx.updateTime(params_.tstep);
  Log::trace() << "IdentityModel<MODEL>::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityModel<MODEL>::finalize(State_ & xx) const {
  Log::trace() << "IdentityModel<MODEL>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_IDENTITYMODEL_H_
