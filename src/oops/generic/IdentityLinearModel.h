/*
 * (C) Copyright 2020-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_IDENTITYLINEARMODEL_H_
#define OOPS_GENERIC_IDENTITYLINEARMODEL_H_

#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/LinearModelBase.h"
#include "oops/interface/Increment.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

class IdentityLinearModelParameters : public LinearModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(IdentityLinearModelParameters, LinearModelParametersBase)

 public:
  oops::RequiredParameter<util::Duration> tstep{"tstep", this};
  oops::RequiredParameter<Variables> vars{"increment variables", this};
  // This does not really belong here but adding as a quick fix
  oops::OptionalParameter<std::string> variableChange{"variable change", this};
};

/// Generic implementation of identity linear model
template <typename MODEL>
class IdentityLinearModel : public LinearModelBase<MODEL> {
  typedef Geometry<MODEL>          Geometry_;
  typedef Increment<MODEL>         Increment_;
  typedef ModelAuxControl<MODEL>   ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL> ModelAuxInc_;
  typedef State<MODEL>             State_;

 public:
  typedef IdentityLinearModelParameters           Parameters_;

  static const std::string classname() {return "oops::IdentityLinearModel";}

  IdentityLinearModel(const Geometry_ &, const IdentityLinearModelParameters &);

/// initialize tangent linear forecast
  void initializeTL(Increment_ &) const override;
/// one tangent linear forecast step
  void stepTL(Increment_ &, const ModelAuxInc_ &) const override;
/// finalize tangent linear forecast
  void finalizeTL(Increment_ &) const override;

/// initialize adjoint forecast
  void initializeAD(Increment_ &) const override;
/// one adjoint forecast step
  void stepAD(Increment_ &, ModelAuxInc_ &) const override;
/// finalize adjoint forecast
  void finalizeAD(Increment_ &) const override;

/// set trajectory
  void setTrajectory(const State_ &, State_ &, const ModelAuxCtl_ &) override;

/// linear model time step
  const util::Duration & timeResolution() const override {return params_.tstep;}
/// linear model variables
  const oops::Variables & variables() const override {return params_.vars;}

 private:
  void print(std::ostream &) const override {}
  const IdentityLinearModelParameters params_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
IdentityLinearModel<MODEL>::IdentityLinearModel(const Geometry_ & resol,
                                                const IdentityLinearModelParameters & params)
  : params_(params) {
  Log::trace() << "IdentityLinearModel<MODEL>::IdentityLinearModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityLinearModel<MODEL>::initializeTL(Increment_ & dx) const {
  Log::trace() << "IdentityLinearModel<MODEL>::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityLinearModel<MODEL>::stepTL(Increment_ & dx, const ModelAuxInc_ & merr) const {
  Log::trace() << "IdentityLinearModel<MODEL>:stepTL Starting " << std::endl;
  dx.updateTime(params_.tstep);
  Log::trace() << "IdentityLinearModel<MODEL>::stepTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityLinearModel<MODEL>::finalizeTL(Increment_ & dx) const {
  Log::trace() << "IdentityLinearModel<MODEL>::finalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityLinearModel<MODEL>::initializeAD(Increment_ & dx) const {
  Log::trace() << "IdentityLinearModel<MODEL>::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityLinearModel<MODEL>::stepAD(Increment_ & dx, ModelAuxInc_ & merr) const {
  Log::trace() << "IdentityLinearModel<MODEL>:stepAD Starting " << std::endl;
  const util::Duration tstep = params_.tstep;
  dx.updateTime(-tstep);
  Log::trace() << "IdentityLinearModel<MODEL>::stepAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityLinearModel<MODEL>::finalizeAD(Increment_ & dx) const {
  Log::trace() << "IdentityLinearModel<MODEL>::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IdentityLinearModel<MODEL>::setTrajectory(const State_ & x, State_ & xlr,
                                               const ModelAuxCtl_ & maux) {
  Log::trace() << "IdentityLinearModel<MODEL>::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_IDENTITYLINEARMODEL_H_
