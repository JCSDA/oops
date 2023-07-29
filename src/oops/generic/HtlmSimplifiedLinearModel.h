/*
 * (C) Copyright 2023 UCAR.
 * (C) Crown copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HTLMSIMPLIFIEDLINEARMODEL_H_
#define OOPS_GENERIC_HTLMSIMPLIFIEDLINEARMODEL_H_

#include "oops/base/LinearModel.h"

namespace oops {

template <typename MODEL>
class HtlmSimplifiedLinearModelParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HtlmSimplifiedLinearModelParameters, Parameters)

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> linearModel{"linear model", this};
};

//------------------------------------------------------------------------------

template <typename MODEL>
class HtlmSimplifiedLinearModel {
  typedef Geometry<MODEL>                               Geometry_;
  typedef HtlmSimplifiedLinearModelParameters<MODEL>    Parameters_;
  typedef Increment<MODEL>                              Increment_;
  typedef LinearModel<MODEL>                            LinearModel_;
  typedef ModelAuxControl<MODEL>                        ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL>                      ModelAuxInc_;
  typedef State<MODEL>                                  State_;

 public:
  HtlmSimplifiedLinearModel(const Parameters_ & params, const Geometry_ & updateGeometry)
  : linearModel_(updateGeometry, params.linearModel) {}

  void forecastSimplifiedTL(Increment_ & dx, const ModelAuxInc_ & merr,
                            const util::Duration & updateTimestep) const {
    linearModel_.forecastTL(dx, merr, updateTimestep);
  }
  void forecastSimplifiedAD(Increment_ & dx, ModelAuxInc_ & merr,
                            const util::Duration & updateTimestep) const {
    linearModel_.forecastAD(dx, merr, updateTimestep);
  }

  void setSimplifiedTrajectory(const State_ & x, State_ & xlr, const ModelAuxCtl_ & maux) {
    linearModel_.setTrajectory(x, xlr, maux);
  }

  const util::Duration & timeResolution() const {return linearModel_.timeResolution();}
  const oops::Variables & variables() const {return linearModel_.variables();}

 private:
  LinearModel_ linearModel_;
};

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMSIMPLIFIEDLINEARMODEL_H_

