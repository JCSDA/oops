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
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "oops/base/LinearObsOperators.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsFilters.h"
#include "oops/base/ObsOperators.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/InterpolatorTraj.h"
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
class Observer : public PostBase<STATE>,
                 public util::Printable {
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef InterpolatorTraj<MODEL>    InterpolatorTraj_;
  typedef LinearObsOperators<MODEL>  LinearObsOperators_;
  typedef ObsAuxControls<MODEL>      ObsAuxCtrls_;
  typedef ObsFilters<MODEL>          ObsFilters_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsOperators<MODEL>        ObsOperators_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef std::vector<boost::shared_ptr<InterpolatorTraj_> > type_vspit;
  typedef boost::shared_ptr<ObsFilters_> PtrFilters_;

 public:
  Observer(const ObsSpaces_ &, const ObsOperators_ &, const ObsAuxCtrls_ &,
           const std::vector<PtrFilters_> filters = std::vector<PtrFilters_>(0),
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
  void print(std::ostream &) const override;

// Obs operator
  const ObsOperators_ & hop_;

// Data
  std::auto_ptr<Observations_> yobs_;
  const ObsAuxCtrls_ & ybias_;

  util::DateTime winbgn_;   //!< Begining of assimilation window
  util::DateTime winend_;   //!< End of assimilation window
  util::DateTime bgn_;      //!< Begining of currently active observations
  util::DateTime end_;      //!< End of currently active observations
  util::Duration hslot_;    //!< Half time slot
  const bool subwindows_;

  std::vector<PtrFilters_> filters_;
  std::vector<Variables> geovars_;  // Variables needed from model (through geovals)
  std::vector<boost::shared_ptr<GeoVaLs_> > gvals_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
Observer<MODEL, STATE>::Observer(const ObsSpaces_ & obsdb,
                                 const ObsOperators_ & hop,
                                 const ObsAuxCtrls_ & ybias,
                                 const std::vector<PtrFilters_> filters,
                                 const util::Duration & tslot, const bool swin)
  : PostBase<STATE>(), hop_(hop),
    yobs_(new Observations_(obsdb)), ybias_(ybias),
    winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()),
    bgn_(winbgn_), end_(winend_), hslot_(tslot/2), subwindows_(swin),
    filters_(filters), geovars_(hop_.size()), gvals_(0)
{
  Log::trace() << "Observer::Observer starting" << std::endl;
  if (filters_.empty()) {
    for (size_t jj = 0; jj < obsdb.size(); ++jj) {
      boost::shared_ptr<ObsFilters_> tmp(new ObsFilters_());
      filters_.push_back(tmp);
    }
  }
  ASSERT(filters_.size() == hop_.size() && ybias_.size() == hop_.size());

  for (size_t jj = 0; jj < hop_.size(); ++jj) {
    geovars_[jj] += hop_[jj].variables();
    geovars_[jj] += ybias_[jj].variables();
    geovars_[jj] += filters_.at(jj)->requiredGeoVaLs();
  }

  Log::trace() << "Observer::Observer done" << std::endl;
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
    boost::shared_ptr<GeoVaLs_> tmp(new GeoVaLs_(hop_[jj].locations(bgn_, end_), geovars_[jj]));
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
    xx.getValues(hop_[jj].locations(t1, t2), geovars_[jj], *gvals_.at(jj));
  }
  Log::trace() << "Observer::doProcessing done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::processTraj(const STATE & xx, type_vspit & traj) const {
  Log::trace() << "Observer::processTraj start" << std::endl;
  util::DateTime t1(xx.validTime()-hslot_);
  util::DateTime t2(xx.validTime()+hslot_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

// Get state variables at obs locations and trajectory
  for (size_t jj = 0; jj < hop_.size(); ++jj) {
    xx.getValues(hop_[jj].locations(t1, t2), geovars_[jj], *gvals_.at(jj),
                 *traj.at(jj));
  }
  Log::trace() << "Observer::processTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::finalizeTraj(const STATE & xx, LinearObsOperators_ & htlad) {
  Log::trace() << "Observer::finalizeTraj start" << std::endl;
  for (size_t jj = 0; jj < htlad.size(); ++jj) {
    htlad[jj].setTrajectory(*gvals_.at(jj), ybias_[jj]);
  }
  this->doFinalize(xx);
  Log::trace() << "Observer::finalizeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::doFinalize(const STATE &) {
  Log::trace() << "Observer::doFinalize start" << std::endl;
  for (size_t jj = 0; jj < hop_.size(); ++jj) {
    filters_.at(jj)->priorFilter(*gvals_.at(jj));
    hop_[jj].simulateObs(*gvals_.at(jj), (*yobs_)[jj], ybias_[jj]);
    filters_.at(jj)->postFilter((*yobs_)[jj]);
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
