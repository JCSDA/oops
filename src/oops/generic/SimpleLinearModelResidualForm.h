/*
 * (C) Copyright 2023 UCAR.
 * (C) Crown copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_SIMPLELINEARMODELRESIDUALFORM_H_
#define OOPS_GENERIC_SIMPLELINEARMODELRESIDUALFORM_H_

#include "oops/generic/SimpleLinearModelMultiresolution.h"

namespace oops {

/// \brief Wrapper for the MODEL-specific LinearModel used as part of the HybridLinearModel, for the
/// case where the resolution of the HybridLinearModel coefficients is not equal to that of the
/// MODEL-specific LinearModel, and the residual form algorithm has been requested.
template <typename MODEL>
class SimpleLinearModelResidualForm : public SimpleLinearModelMultiresolution<MODEL> {
  typedef Geometry<MODEL>             Geometry_;
  typedef Increment<MODEL>            Increment_;
  typedef LinearModel<MODEL>          LinearModel_;
  typedef ModelAuxControl<MODEL>      ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL>    ModelAuxInc_;
  typedef State<MODEL>                State_;

 public:
  SimpleLinearModelResidualForm(const eckit::Configuration &, const Geometry_ &,
                                const Variables &, const util::DateTime &);
  void forecastTL(Increment_ &, const ModelAuxInc_ &, const util::Duration &) const override;
  void forecastAD(Increment_ &, ModelAuxInc_ &, const util::Duration &) const override;

 private:
  mutable Increment_ dxLowResCopy_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
SimpleLinearModelResidualForm<MODEL>::SimpleLinearModelResidualForm(
                                                                const eckit::Configuration & config,
                                                                const Geometry_ & updateGeometry,
                                                                const Variables & vars,
                                                                const util::DateTime & windowBegin)
: SimpleLinearModelMultiresolution<MODEL>(config, updateGeometry),
  dxLowResCopy_(this->geometry_, vars, windowBegin) {}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SimpleLinearModelResidualForm<MODEL>::forecastTL(Increment_ & dx,
                                                      const ModelAuxInc_ & merr,
                                                      const util::Duration & updateTstep) const {
  Increment_ dxLowRes(this->geometry_, dx);
  dxLowResCopy_ = dxLowRes;
  this->linearModel_->forecastTL(dxLowRes, merr, updateTstep);
  dxLowResCopy_.updateTime(updateTstep);
  dxLowRes -= dxLowResCopy_;
  dx.updateTime(updateTstep);
  dx += Increment_(this->updateGeometry_, dxLowRes);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SimpleLinearModelResidualForm<MODEL>::forecastAD(Increment_ & dx,
                                                      ModelAuxInc_ & merr,
                                                      const util::Duration & updateTstep) const {
  Increment_ dxLowRes(this->geometry_, dx, true);
  dxLowResCopy_ = dxLowRes;
  this->linearModel_->forecastAD(dxLowRes, merr, updateTstep);
  dxLowResCopy_.updateTime(-updateTstep);
  dxLowRes -= dxLowResCopy_;
  dx.updateTime(-updateTstep);
  dx += Increment_(this->updateGeometry_, dxLowRes, true);
}

}  // namespace oops

#endif  // OOPS_GENERIC_SIMPLELINEARMODELRESIDUALFORM_H_
