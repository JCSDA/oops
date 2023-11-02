/*
 * (C) Copyright 2023 UCAR.
 * (C) Crown copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HTLMSIMPLIFIEDLINEARMODEL_H_
#define OOPS_GENERIC_HTLMSIMPLIFIEDLINEARMODEL_H_

#include <memory>

#include "oops/base/LinearModel.h"
#include "oops/base/TrajectorySaver.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class HtlmSimplifiedLinearModel {
  typedef Geometry<MODEL>                               Geometry_;
  typedef Increment<MODEL>                              Increment_;
  typedef LinearModel<MODEL>                            LinearModel_;
  typedef ModelAuxControl<MODEL>                        ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL>                      ModelAuxInc_;
  typedef State<MODEL>                                  State_;

 public:
  HtlmSimplifiedLinearModel(const eckit::Configuration & config, const Geometry_ & updateGeometry)
  : config_(eckit::LocalConfiguration(config, "linear model")), updateGeometry_(updateGeometry) {
    linearModel_ = std::make_shared<LinearModel_>(updateGeometry, config_);
  }

  void forecastSimplifiedTL(Increment_ & dx, const ModelAuxInc_ & merr,
                            const util::Duration & updateTimestep) const {
    linearModel_->forecastTL(dx, merr, updateTimestep);
  }
  void forecastSimplifiedAD(Increment_ & dx, ModelAuxInc_ & merr,
                            const util::Duration & updateTimestep) const {
    linearModel_->forecastAD(dx, merr, updateTimestep);
  }

  void setUpTrajectorySaver(PostProcessor<State_> & pp, const ModelAuxCtl_ & maux) {
    PostProcessorTLAD<MODEL> pptraj;
    pp.enrollProcessor(new TrajectorySaver<MODEL>(config_, updateGeometry_, maux, linearModel_,
                                                  pptraj));
  }

  void setSimplifiedTrajectory(const State_ & x, State_ & xlr, const ModelAuxCtl_ & maux) {
    linearModel_->setTrajectory(x, xlr, maux);
  }

  const util::Duration & timeResolution() const {return linearModel_->timeResolution();}
  const oops::Variables & variables() const {return linearModel_->variables();}

 private:
  const eckit::LocalConfiguration config_;
  const Geometry_ & updateGeometry_;
  std::shared_ptr<LinearModel_> linearModel_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMSIMPLIFIEDLINEARMODEL_H_

