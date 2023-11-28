/*
 * (C) Copyright 2023 UCAR.
 * (C) Crown copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_SIMPLELINEARMODEL_H_
#define OOPS_GENERIC_SIMPLELINEARMODEL_H_

#include <memory>

#include "oops/base/LinearModel.h"
#include "oops/base/TrajectorySaver.h"

namespace oops {

/// \brief Wrapper for the MODEL-specific LinearModel used as part of the HybridLinearModel. This
/// class is for the case where the resolution of the HybridLinearModel coefficients is equal to
/// that of the MODEL-specific LinearModel; derived classes handle the converse.
template <typename MODEL>
class SimpleLinearModel {
  typedef Geometry<MODEL>             Geometry_;
  typedef Increment<MODEL>            Increment_;
  typedef LinearModel<MODEL>          LinearModel_;
  typedef ModelAuxControl<MODEL>      ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL>    ModelAuxInc_;
  typedef PostProcessorTLAD<MODEL>    PostProcessorTLAD_;
  typedef State<MODEL>                State_;
  typedef TrajectorySaver<MODEL>      TrajectorySaver_;

 public:
  SimpleLinearModel(const eckit::Configuration & config, const Geometry_ & updateGeometry);
  virtual void forecastTL(Increment_ &, const ModelAuxInc_ &, const util::Duration &) const;
  virtual void forecastAD(Increment_ &, ModelAuxInc_ &, const util::Duration &) const;
  virtual void setUpTrajectorySaver(PostProcessor<State_> &, const ModelAuxCtl_ &);
  virtual void setTrajectory(const State_ &, State_ &, const ModelAuxCtl_ &);
  const util::Duration & timeResolution() const {return linearModel_->timeResolution();}
  const oops::Variables & variables() const {return linearModel_->variables();}

 protected:
  const eckit::LocalConfiguration lmConfig_;
  const Geometry_ & updateGeometry_;
  std::shared_ptr<LinearModel_> linearModel_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
SimpleLinearModel<MODEL>::SimpleLinearModel(const eckit::Configuration & config,
                                            const Geometry_ & updateGeometry)
: lmConfig_(config, "linear model"), updateGeometry_(updateGeometry) {
  if (!config.has("geometry")) {
    linearModel_ = std::make_shared<LinearModel_>(updateGeometry, lmConfig_);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SimpleLinearModel<MODEL>::forecastTL(Increment_ & dx,
                                          const ModelAuxInc_ & merr,
                                          const util::Duration & updateTstep) const {
  linearModel_->forecastTL(dx, merr, updateTstep);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SimpleLinearModel<MODEL>::forecastAD(Increment_ & dx,
                                          ModelAuxInc_ & merr,
                                          const util::Duration & updateTstep) const {
  linearModel_->forecastAD(dx, merr, updateTstep);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SimpleLinearModel<MODEL>::setUpTrajectorySaver(PostProcessor<State_> & pp,
                                                    const ModelAuxCtl_ & maux) {
  PostProcessorTLAD_ pptraj;
  pp.enrollProcessor(
    std::make_shared<TrajectorySaver_>(lmConfig_, updateGeometry_, maux, linearModel_, pptraj));
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SimpleLinearModel<MODEL>::setTrajectory(const State_ & x, State_ & xlr,
                                             const ModelAuxCtl_ & maux) {
  linearModel_->setTrajectory(x, xlr, maux);
}

}  // namespace oops

#endif  // OOPS_GENERIC_SIMPLELINEARMODEL_H_
