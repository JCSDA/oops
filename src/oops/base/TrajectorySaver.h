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

#include <memory>

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

// -----------------------------------------------------------------------------

template <typename MODEL>
class TrajectorySaver : public PostBase<State<MODEL> > {
  typedef Geometry<MODEL>          Geometry_;
  typedef LinearModel<MODEL>       LinearModel_;
  typedef ModelAuxControl<MODEL>   ModelAux_;
  typedef PostProcessorTLAD<MODEL> PPTLAD_;
  typedef State<MODEL>             State_;

 public:
  TrajectorySaver(const eckit::Configuration &, const Geometry_ &,
                  const ModelAux_ &, boost::ptr_vector<LinearModel_> &, PPTLAD_);
  TrajectorySaver(const eckit::Configuration &, const Geometry_ &, PPTLAD_);
  ~TrajectorySaver() {}

 private:
  const Geometry_    resol_;
  PPTLAD_ pptraj_;
  const bool fourd_;
  const eckit::LocalConfiguration   tlConf_;
  std::unique_ptr<const ModelAux_>  lrBias_;
  boost::ptr_vector<LinearModel_> * tlm_;
  LinearModel_ *     subtlm_;

  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const State_ &) override;
  void doFinalize(const State_ &) override;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
TrajectorySaver<MODEL>::TrajectorySaver(const eckit::Configuration & conf,
                                        const Geometry_ & resol,
                                        const ModelAux_ & bias,
                                        boost::ptr_vector<LinearModel_> & tlm,
                                        PPTLAD_ pptraj):
  PostBase<State_>(conf),
  resol_(resol), pptraj_(pptraj), fourd_(true),
  tlConf_(conf), lrBias_(new ModelAux_(resol, bias)), tlm_(&tlm), subtlm_(0)
{
  Log::trace() << "TrajectorySaver::TrajectorySaver 4D" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
TrajectorySaver<MODEL>::TrajectorySaver(const eckit::Configuration & conf,
                                        const Geometry_ & resol, PPTLAD_ pptraj):
  PostBase<State_>(conf),
  resol_(resol), pptraj_(pptraj), fourd_(false),
  tlConf_(), lrBias_(), tlm_(nullptr), subtlm_(0)
{
  Log::trace() << "TrajectorySaver::TrajectorySaver 3D" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void TrajectorySaver<MODEL>::doInitialize(const State_ & x0,
                                          const util::DateTime & end,
                                          const util::Duration & step) {
  Log::trace() << "TrajectorySaver::doInitialize start" << std::endl;
  if (fourd_) subtlm_ = new LinearModel_(resol_, tlConf_);
  State_ xlr(resol_, x0);
  pptraj_.initializeTraj(xlr, end, step);
  Log::trace() << "TrajectorySaver::doInitialize done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void TrajectorySaver<MODEL>::doProcessing(const State_ & xx) {
  Log::trace() << "TrajectorySaver::doProcessing start" << std::endl;
  State_ xlr(resol_, xx);
  if (fourd_) {
    ASSERT(subtlm_ != 0);
    subtlm_->setTrajectory(xx, xlr, *lrBias_);
  }
  pptraj_.processTraj(xlr);
  Log::trace() << "TrajectorySaver::doProcessing done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void TrajectorySaver<MODEL>::doFinalize(const State_ & xx) {
  Log::trace() << "TrajectorySaver::doFinalize start" << std::endl;
  State_ xlr(resol_, xx);
  if (fourd_) tlm_->push_back(subtlm_);
  pptraj_.finalizeTraj(xlr);
  subtlm_ = 0;
  Log::trace() << "TrajectorySaver::doFinalize done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_TRAJECTORYSAVER_H_
