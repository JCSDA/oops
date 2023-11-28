/*
 * (C) Copyright 2023 UCAR.
 * (C) Crown copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_SIMPLELINEARMODELMULTIRESOLUTION_H_
#define OOPS_GENERIC_SIMPLELINEARMODELMULTIRESOLUTION_H_

#include <memory>

#include "oops/generic/SimpleLinearModel.h"

namespace oops {

/// \brief Wrapper for the MODEL-specific LinearModel used as part of the HybridLinearModel, for the
/// case where the resolution of the HybridLinearModel coefficients is not equal to that of the
/// MODEL-specific LinearModel.
template <typename MODEL>
class SimpleLinearModelMultiresolution : public SimpleLinearModel<MODEL> {
  typedef Geometry<MODEL>             Geometry_;
  typedef Increment<MODEL>            Increment_;
  typedef LinearModel<MODEL>          LinearModel_;
  typedef ModelAuxControl<MODEL>      ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL>    ModelAuxInc_;
  typedef PostProcessorTLAD<MODEL>    PostProcessorTLAD_;
  typedef SimpleLinearModel<MODEL>    SimpleLinearModel_;
  typedef State<MODEL>                State_;
  typedef TrajectorySaver<MODEL>      TrajectorySaver_;

 public:
  SimpleLinearModelMultiresolution(const eckit::Configuration &, const Geometry_ &);
  void forecastTL(Increment_ &, const ModelAuxInc_ &, const util::Duration &) const override;
  void forecastAD(Increment_ &, ModelAuxInc_ &, const util::Duration &) const override;
  void setUpTrajectorySaver(PostProcessor<State_> &, const ModelAuxCtl_ &) final;
  void setTrajectory(const State_ &, State_ &, const ModelAuxCtl_ &) final;

 protected:
  const Geometry_ geometry_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
SimpleLinearModelMultiresolution<MODEL>::SimpleLinearModelMultiresolution(
                                                                const eckit::Configuration & config,
                                                                const Geometry_ & updateGeometry)
: SimpleLinearModel_(config, updateGeometry),
  geometry_(eckit::LocalConfiguration(config, "geometry"), updateGeometry.getComm()) {
  this->linearModel_ = std::make_shared<LinearModel_>(geometry_, this->lmConfig_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SimpleLinearModelMultiresolution<MODEL>::forecastTL(Increment_ & dx,
                                                         const ModelAuxInc_ & merr,
                                                         const util::Duration & updateTstep) const {
  Increment_ dxLowRes(geometry_, dx);
  this->linearModel_->forecastTL(dxLowRes, merr, updateTstep);
  dx = Increment(this->updateGeometry_, dxLowRes);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SimpleLinearModelMultiresolution<MODEL>::forecastAD(Increment_ & dx,
                                                         ModelAuxInc_ & merr,
                                                         const util::Duration & updateTstep) const {
  Increment_ dxLowRes(geometry_, dx, true);
  this->linearModel_->forecastAD(dxLowRes, merr, updateTstep);
  dx = Increment(this->updateGeometry_, dxLowRes, true);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SimpleLinearModelMultiresolution<MODEL>::setUpTrajectorySaver(PostProcessor<State_> & pp,
                                                                   const ModelAuxCtl_ & maux) {
  PostProcessorTLAD_ pptraj;
  pp.enrollProcessor(std::make_shared<TrajectorySaver_>(
    this->lmConfig_, geometry_, maux, this->linearModel_, pptraj));
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SimpleLinearModelMultiresolution<MODEL>::setTrajectory(const State_ & x,
                                                            State_ & xlr,
                                                            const ModelAuxCtl_ & maux) {
  State_ xLowRes(geometry_, x);
  this->linearModel_->setTrajectory(x, xLowRes, maux);
}

}  // namespace oops

#endif  // OOPS_GENERIC_SIMPLELINEARMODELMULTIRESOLUTION_H_
