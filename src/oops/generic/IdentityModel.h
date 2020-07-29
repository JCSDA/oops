/*
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_IDENTITYMODEL_H_
#define OOPS_GENERIC_IDENTITYMODEL_H_

#include <string>

#include "oops/base/ModelBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// Generic implementation of identity model
template <typename MODEL>
class IdentityModel : public ModelBase<MODEL> {
  typedef typename MODEL::Geometry          Geometry_;
  typedef typename MODEL::ModelAuxControl   ModelAux_;
  typedef typename MODEL::State             State_;

 public:
  static const std::string classname() {return "oops::IdentityModel";}

  IdentityModel(const Geometry_ &, const eckit::Configuration &);

/// initialize forecast
  void initialize(State_ &) const override;
/// one forecast step
  void step(State_ &, const ModelAux_ &) const override;
/// finalize forecast
  void finalize(State_ &) const override;

/// model time step
  const util::Duration & timeResolution() const override {return tstep_;}
/// model variables
  const oops::Variables & variables() const override {return vars_;}

 private:
  void print(std::ostream &) const override {}
  const util::Duration tstep_;
  const oops::Variables vars_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
IdentityModel<MODEL>::IdentityModel(const Geometry_ & resol, const eckit::Configuration & conf)
  : tstep_(util::Duration(conf.getString("tstep"))), vars_(conf, "state variables") {
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
  xx.updateTime(tstep_);
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
