/*
 * (C) Copyright 2020-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/ModelBase.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// Generic implementation of identity model
template <typename MODEL>
class IdentityModel : public ModelBase<MODEL> {
  typedef Geometry<MODEL>          Geometry_;
  typedef ModelAuxControl<MODEL>   ModelAux_;
  typedef State<MODEL>             State_;

 public:
  static const std::string classname() {return "oops::IdentityModel";}

  IdentityModel(const Geometry_ &, const eckit::Configuration &);

/// initialize forecast
  void initialize(State_ &) const override {}
/// one forecast step
  void step(State_ &, const ModelAux_ &) const override;
/// finalize forecast
  void finalize(State_ &) const override {}

/// model time step
  const util::Duration & timeResolution() const override {return tstep_;}

 private:
  void print(std::ostream &) const override {}
  const util::Duration tstep_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
IdentityModel<MODEL>::IdentityModel(const Geometry_ & resol, const eckit::Configuration & config)
  : tstep_(config.getString("tstep")) {
  Log::trace() << "IdentityModel<MODEL>::IdentityModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityModel<MODEL>::step(State_ & xx, const ModelAux_ & merr) const {
  Log::trace() << "IdentityModel<MODEL>:step Starting " << std::endl;
  xx.updateTime(tstep_);
  Log::trace() << "IdentityModel<MODEL>::step done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
