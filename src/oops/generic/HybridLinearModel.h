/*
 * (C) Copyright 2022-2023 UCAR.
 * (C) Crown copyright 2022-2023 Met Office.
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
#include "oops/generic/HybridLinearModelCoeffs.h"
#include "oops/generic/LinearModelBase.h"
#include "oops/generic/SimpleLinearModelResidualForm.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

//------------------------------------------------------------------------------

template <typename MODEL>
class HybridLinearModel : public LinearModelBase<MODEL> {
  typedef Geometry<MODEL>                         Geometry_;
  typedef HtlmCalculator<MODEL>                   HtlmCalculator_;
  typedef HtlmEnsemble<MODEL>                     HtlmEnsemble_;
  typedef HybridLinearModelCoeffs<MODEL>          HybridLinearModelCoeffs_;
  typedef Increment<MODEL>                        Increment_;
  typedef ModelAuxControl<MODEL>                  ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL>                ModelAuxInc_;
  typedef SimpleLinearModel<MODEL>                SimpleLinearModel_;
  typedef SimpleLinearModelMultiresolution<MODEL> SLMMultiresolution_;
  typedef SimpleLinearModelResidualForm<MODEL>    SLMResidualForm_;
  typedef State<MODEL>                            State_;

 public:
  static const std::string classname() {return "oops::HybridLinearModel";}

  HybridLinearModel(const Geometry_ &, const eckit::Configuration &);

  void initializeTL(Increment_ &) const override {}
  void stepTL(Increment_ &, const ModelAuxInc_ &) const override;
  void finalizeTL(Increment_ &) const override {}

  void initializeAD(Increment_ &) const override {}
  void stepAD(Increment_ &, ModelAuxInc_ &) const override;
  void finalizeAD(Increment_ &) const override {}

  void setTrajectory(const State_ &, State_ &, const ModelAuxCtl_ &) override;

  const util::Duration & timeResolution() const override {return updateTstep_;}

 private:
  void print(std::ostream &) const override {}
  const util::Duration updateTstep_;
  const Variables vars_;  // superset of coeffs_::updateVars_
  std::unique_ptr<SimpleLinearModel_> simpleLinearModel_;
  HybridLinearModelCoeffs_ coeffs_;
};

// ------------------------------------------------------------------------------

template<typename MODEL>
HybridLinearModel<MODEL>::HybridLinearModel(const Geometry_ & updateGeometry,
                                            const eckit::Configuration & config)
: updateTstep_(config.getString("update tstep")), vars_(config, "variables"),
  coeffs_(config.getSubConfiguration("coefficients"), updateGeometry, updateTstep_) {
  Log::trace() << "HybridLinearModel<MODEL>::HybridLinearModel starting" << std::endl;
  // Set up simpleLinearModel_
  const eckit::LocalConfiguration slmConf(config, "simple linear model");
  const util::DateTime wBgn(config.getSubConfiguration("coefficients").getString("window begin"));
  if (slmConf.getBool("residual form", false)) {
    simpleLinearModel_ = std::make_unique<SLMResidualForm_>(slmConf, updateGeometry, vars_, wBgn);
  } else if (slmConf.has("geometry")) {
    simpleLinearModel_ = std::make_unique<SLMMultiresolution_>(slmConf, updateGeometry);
  } else {
    simpleLinearModel_ = std::make_unique<SimpleLinearModel_>(slmConf, updateGeometry);
  }
  if (updateTstep_ % simpleLinearModel_->timeResolution() != 0) {
    ABORT("HybridLinearModel<MODEL>::HybridLinearModel: "
          "update tstep is not a multiple of simple linear model tstep");
  }
  // Obtain coefficients from file or by generating them
  coeffs_.obtain(*simpleLinearModel_, vars_);
  Log::trace() << "HybridLinearModel<MODEL>::HybridLinearModel done" << std::endl;
}

// ------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModel<MODEL>::stepTL(Increment_ & dx, const ModelAuxInc_ & merr) const {
  Log::trace() << "HybridLinearModel<MODEL>::stepTL() starting" << std::endl;
  simpleLinearModel_->forecastTL(dx, merr, updateTstep_);
  coeffs_.updateIncTL(dx);
  Log::trace() << "HybridLinearModel<MODEL>::stepTL() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModel<MODEL>::stepAD(Increment_ & dx, ModelAuxInc_ & merr) const {
  Log::trace() << "HybridLinearModel<MODEL>::stepAD() starting" << std::endl;
  coeffs_.updateIncAD(dx);
  simpleLinearModel_->forecastAD(dx, merr, updateTstep_);
  Log::trace() << "HybridLinearModel<MODEL>::stepAD() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModel<MODEL>::setTrajectory(const State_ & x, State_ & xlr,
                                             const ModelAuxCtl_ & maux) {
  Log::trace() << "HybridLinearModel<MODEL>::setTrajectory() starting" << std::endl;
  simpleLinearModel_->setTrajectory(x , xlr, maux);
  Log::trace() << "HybridLinearModel<MODEL>::setTrajectory() done" << std::endl;
}

//------------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_HYBRIDLINEARMODEL_H_
