/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_TRAJECTORYSAVER_H_
#define OOPS_BASE_TRAJECTORYSAVER_H_

#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/PostBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/LinearModel.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

/// Save trajectory during forecast run.

template <typename MODEL>
class TrajectorySaver : public PostBase<State<MODEL> > {
  typedef Geometry<MODEL>        Geometry_;
  typedef LinearModel<MODEL>     LinearModel_;
  typedef ModelAuxControl<MODEL> ModelAux_;
  typedef State<MODEL>           State_;

 public:
  TrajectorySaver(const State_ &, const eckit::Configuration &,
                  const Geometry_ &, const ModelAux_ &,
                  boost::ptr_vector<LinearModel_> &);
  ~TrajectorySaver() {}

 private:
  const Geometry_    resol_;
  const eckit::LocalConfiguration tlConf_;
  const ModelAux_    lrBias_;
  boost::ptr_vector<LinearModel_> & tlm_;
  LinearModel_ *     subtlm_;
  State_             xlr_;

  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const State_ &) override;
  void doFinalize(const State_ &) override;
};

// ====================================================================================

template <typename MODEL>
TrajectorySaver<MODEL>::TrajectorySaver(const State_ & xx,
                                        const eckit::Configuration & conf,
                                        const Geometry_ & resol,
                                        const ModelAux_ & bias,
                                        boost::ptr_vector<LinearModel_> & tlm):
  PostBase<State_>(conf),
  resol_(resol), tlConf_(conf), lrBias_(resol, bias),
  tlm_(tlm), subtlm_(0), xlr_(resol, xx)
{}
// -----------------------------------------------------------------------------
template <typename MODEL>
void TrajectorySaver<MODEL>::doInitialize(const State_ &,
                                          const util::DateTime &, const util::Duration &) {
  subtlm_ = new LinearModel_(resol_, tlConf_);
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void TrajectorySaver<MODEL>::doProcessing(const State_ & xx) {
  ASSERT(subtlm_ != 0);
  subtlm_->setTrajectory(xx, xlr_, lrBias_);
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void TrajectorySaver<MODEL>::doFinalize(const State_ &) {
  tlm_.push_back(subtlm_);
  subtlm_ = 0;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_TRAJECTORYSAVER_H_
