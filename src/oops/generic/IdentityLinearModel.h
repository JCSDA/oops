/*
 * (C) Copyright 2020-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/LinearModelBase.h"
#include "oops/interface/Increment.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Generic implementation of identity linear model
template <typename MODEL>
class IdentityLinearModel : public LinearModelBase<MODEL> {
  typedef Geometry<MODEL>          Geometry_;
  typedef Increment<MODEL>         Increment_;
  typedef ModelAuxControl<MODEL>   ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL> ModelAuxInc_;
  typedef State<MODEL>             State_;

 public:
  static const std::string classname() {return "oops::IdentityLinearModel";}

  IdentityLinearModel(const Geometry_ &, const eckit::Configuration &);

/// initialize tangent linear forecast
  void initializeTL(Increment_ &) const override {}
/// one tangent linear forecast step
  void stepTL(Increment_ &, const ModelAuxInc_ &) const override;
/// finalize tangent linear forecast
  void finalizeTL(Increment_ &) const override {}

/// initialize adjoint forecast
  void initializeAD(Increment_ &) const override {}
/// one adjoint forecast step
  void stepAD(Increment_ &, ModelAuxInc_ &) const override;
/// finalize adjoint forecast
  void finalizeAD(Increment_ &) const override {}

/// set trajectory
  void setTrajectory(const State_ &, State_ &, const ModelAuxCtl_ &) override {}

/// linear model time step
  const util::Duration & timeResolution() const override {return tstep_;}
/// linear model variables
  const oops::Variables & variables() const override {return vars_;}

 private:
  void print(std::ostream &) const override {}
  const util::Duration tstep_;
  const oops::Variables vars_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
IdentityLinearModel<MODEL>::IdentityLinearModel(const Geometry_ & resol,
                                                const eckit::Configuration & config)
  : tstep_(config.getString("tstep")), vars_(config, "increment variables") {
  Log::trace() << "IdentityLinearModel<MODEL>::IdentityLinearModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityLinearModel<MODEL>::stepTL(Increment_ & dx, const ModelAuxInc_ &) const {
  Log::trace() << "IdentityLinearModel<MODEL>:stepTL Starting " << std::endl;
  dx.updateTime(tstep_);
  Log::trace() << "IdentityLinearModel<MODEL>::stepTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityLinearModel<MODEL>::stepAD(Increment_ & dx, ModelAuxInc_ &) const {
  Log::trace() << "IdentityLinearModel<MODEL>:stepAD Starting " << std::endl;
  dx.updateTime(-tstep_);
  Log::trace() << "IdentityLinearModel<MODEL>::stepAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
