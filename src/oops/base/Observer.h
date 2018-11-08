/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVER_H_
#define OOPS_BASE_OBSERVER_H_

#include <memory>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "oops/base/LinearObsOperators.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsFilters.h"
#include "oops/base/ObsOperators.h"
#include "oops/base/PostBase.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/InterpolatorTraj.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

/// Computes observation equivalent during model run.

// Sub-windows knowledge could be removed if vector of obs was used in
// weak constraint 4D-Var. YT

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
class Observer : public util::Printable, public PostBase<STATE> {
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef InterpolatorTraj<MODEL>    InterpolatorTraj_;
  typedef LinearObsOperators<MODEL>  LinearObsOperators_;
  typedef ObsAuxControl<MODEL>       ObsAuxCtrl_;
  typedef ObsFilters<MODEL>          ObsFilters_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsOperators<MODEL>        ObsOperators_;
  typedef std::vector<boost::shared_ptr<InterpolatorTraj_> > vspit;

 public:
  Observer(const ObsOperators_ &, const ObsAuxCtrl_ &,
           const std::vector<ObsFilters_> &,
           const util::Duration & tslot = util::Duration(0), const bool subwin = false);
  ~Observer() {}

  Observations_ * release() {return yobs_.release();}

  void processTraj(const STATE &, std::vector<boost::shared_ptr<InterpolatorTraj_> > &) const;
  void finalizeTraj(const STATE &, LinearObsOperators_ &);

 private:
// Methods
  void doInitialize(const STATE &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const STATE &) override;
  void doFinalize(const STATE &) override;
  void print(std::ostream &) const;

// Obs operator
  const ObsOperators_ & hop_;

// Data
  std::auto_ptr<Observations_> yobs_;
  const ObsAuxCtrl_ & ybias_;

  util::DateTime winbgn_;   //!< Begining of assimilation window
  util::DateTime winend_;   //!< End of assimilation window
  util::DateTime bgn_;      //!< Begining of currently active observations
  util::DateTime end_;      //!< End of currently active observations
  util::Duration hslot_;    //!< Half time slot
  const bool subwindows_;

  std::vector<boost::shared_ptr<GeoVaLs_> > gvals_;
  std::vector<ObsFilters_> filters_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
Observer<MODEL, STATE>::Observer(const ObsOperators_ & hop,
                                 const ObsAuxCtrl_ & ybias,
                                 const std::vector<ObsFilters_> & filters,
                                 const util::Duration & tslot, const bool swin)
  : PostBase<STATE>(), hop_(hop),
    yobs_(new Observations_(hop_)), ybias_(ybias),
    winbgn_(hop_.windowStart()), winend_(hop_.windowEnd()),
    bgn_(winbgn_), end_(winend_), hslot_(tslot/2), subwindows_(swin),
    gvals_(0), filters_(filters)
{
  Log::trace() << "Observer::Observer" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::doInitialize(const STATE & xx,
                                          const util::DateTime & end,
                                          const util::Duration & tstep) {
  Log::trace() << "Observer::doInitialize start" << std::endl;
  const util::DateTime bgn(xx.validTime());
  if (hslot_ == util::Duration(0)) hslot_ = tstep/2;
  if (subwindows_) {
    if (bgn == end) {
      bgn_ = bgn - hslot_;
      end_ = end + hslot_;
    } else {
      bgn_ = bgn;
      end_ = end;
    }
  }
  if (bgn_ < winbgn_) bgn_ = winbgn_;
  if (end_ > winend_) end_ = winend_;

  for (size_t jj = 0; jj < hop_.size(); ++jj) {
    boost::shared_ptr<GeoVaLs_> tmp(new GeoVaLs_(hop_[jj].locations(bgn_, end_),
                                                 hop_[jj].variables()));
    gvals_.push_back(tmp);
  }
  Log::trace() << "Observer::doInitialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::doProcessing(const STATE & xx) {
  Log::trace() << "Observer::doProcessing start" << std::endl;
  util::DateTime t1(xx.validTime()-hslot_);
  util::DateTime t2(xx.validTime()+hslot_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

// Get state variables at obs locations
  for (size_t jj = 0; jj < hop_.size(); ++jj) {
    xx.getValues(hop_[jj].locations(t1, t2), hop_[jj].variables(), *gvals_.at(jj));
  }
  Log::trace() << "Observer::doProcessing done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::processTraj(const STATE & xx, vspit & traj) const {
  Log::trace() << "Observer::processTraj start" << std::endl;
  util::DateTime t1(xx.validTime()-hslot_);
  util::DateTime t2(xx.validTime()+hslot_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

// Get state variables at obs locations and trajectory
  for (size_t jj = 0; jj < hop_.size(); ++jj) {
    xx.getValues(hop_[jj].locations(t1, t2), hop_[jj].variables(), *gvals_.at(jj),
                 *traj.at(jj));
  }
  Log::trace() << "Observer::processTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::finalizeTraj(const STATE & xx, LinearObsOperators_ & htlad) {
  Log::trace() << "Observer::finalizeTraj start" << std::endl;
  for (size_t jj = 0; jj < htlad.size(); ++jj) {
    htlad[jj].setTrajectory(*gvals_.at(jj), ybias_);
  }
  this->doFinalize(xx);
  Log::trace() << "Observer::finalizeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::doFinalize(const STATE &) {
  Log::trace() << "Observer::doFinalize start" << std::endl;
  for (size_t jj = 0; jj < hop_.size(); ++jj) {
    filters_[jj].priorFilter(*gvals_.at(jj));
    hop_[jj].simulateObs(*gvals_.at(jj), (*yobs_)[jj], ybias_);
    filters_[jj].postFilter((*yobs_)[jj]);
  }
  gvals_.clear();
  Log::trace() << "Observer::doFinalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::print(std::ostream &) const {}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVER_H_
