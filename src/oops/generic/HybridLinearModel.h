/*
 * (C) Copyright 2021-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HYBRIDLINEARMODEL_H_
#define OOPS_GENERIC_HYBRIDLINEARMODEL_H_

#include <memory>
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

template <typename MODEL>
class HybridLinearModelParameters : public LinearModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(HybridLinearModelParameters, LinearModelParametersBase)

 public:
  oops::RequiredParameter<util::Duration> tstep{"tstep", this};
  oops::RequiredParameter<Variables> vars{"increment variables", this};
  // This does not really belong here but adding as a quick fix
  oops::OptionalParameter<std::string> variableChange{"variable change", this};
  // Option to specify which TLM will be used to generate ensemble for coefficient calculation
  oops::RequiredParameter<eckit::LocalConfiguration> simplifiedTLM{"simplified linear model", this};
};

/// Generic implementation of Hybrid linear model
template <typename MODEL>
class HybridLinearModel : public LinearModelBase<MODEL> {
  typedef Geometry<MODEL>          Geometry_;
  typedef Increment<MODEL>         Increment_;
  typedef ModelAuxControl<MODEL>   ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL> ModelAuxInc_;
  typedef State<MODEL>             State_;
  typedef LinearModelBase<MODEL>              LinearModelBase_;
  typedef LinearModelFactory<MODEL>           LinearModelFactory_;
  typedef LinearModelParametersWrapper<MODEL> LinearModelParametersWrapper_;

 public:
  typedef HybridLinearModelParameters<MODEL>           Parameters_;

  static const std::string classname() {return "oops::HybridLinearModel";}

  HybridLinearModel(const Geometry_ &, const HybridLinearModelParameters<MODEL> &);

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
  const HybridLinearModelParameters<MODEL> params_;
  std::unique_ptr<LinearModelBase_> simplifiedLinearModel_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
HybridLinearModel<MODEL>::HybridLinearModel(const Geometry_ & resol,
                                                const HybridLinearModelParameters<MODEL> & params)
  : params_(params) {
    // instantiate the simplified tangent linear model.
    LinearModelParametersWrapper_ simplifiedTLMParameters;
    simplifiedTLMParameters.validateAndDeserialize(params.simplifiedTLM);
    simplifiedLinearModel_.reset(LinearModelFactory_::create(
          resol, simplifiedTLMParameters.linearModelParameters));
  Log::trace() << "HybridLinearModel<MODEL>::HybridLinearModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModel<MODEL>::initializeTL(Increment_ & dx) const {
  Log::trace() << "HybridLinearModel<MODEL>::initializeTL done" << std::endl;
  simplifiedLinearModel_->initializeTL(dx);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModel<MODEL>::stepTL(Increment_ & dx, const ModelAuxInc_ & merr) const {
  Log::trace() << "HybridLinearModel<MODEL>:stepTL Starting " << std::endl;
  simplifiedLinearModel_->stepTL(dx, merr);
  Log::trace() << "HybridLinearModel<MODEL>::stepTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModel<MODEL>::finalizeTL(Increment_ & dx) const {
  Log::trace() << "HybridLinearModel<MODEL>::finalizeTL done" << std::endl;
  simplifiedLinearModel_->finalizeTL(dx);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModel<MODEL>::initializeAD(Increment_ & dx) const {
  Log::trace() << "HybridLinearModel<MODEL>::initializeAD done" << std::endl;
  simplifiedLinearModel_->initializeAD(dx);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModel<MODEL>::stepAD(Increment_ & dx, ModelAuxInc_ & merr) const {
  Log::trace() << "HybridLinearModel<MODEL>:stepAD Starting " << std::endl;
  simplifiedLinearModel_->stepAD(dx , merr);
  Log::trace() << "HybridLinearModel<MODEL>::stepAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModel<MODEL>::finalizeAD(Increment_ & dx) const {
  Log::trace() << "HybridLinearModel<MODEL>::finalizeAD done" << std::endl;
  simplifiedLinearModel_->finalizeAD(dx);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModel<MODEL>::setTrajectory(const State_ & x, State_ & xlr,
                                               const ModelAuxCtl_ & maux) {
  Log::trace() << "HybridLinearModel<MODEL>::finalizeAD done" << std::endl;
  simplifiedLinearModel_->setTrajectory(x , xlr, maux);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_HYBRIDLINEARMODEL_H_
