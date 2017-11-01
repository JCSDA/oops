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

#include <boost/scoped_ptr.hpp>

#include "oops/base/Observations.h"
#include "oops/base/PostBase.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ModelAtLocations.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsOperator.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

/// Computes observation equivalent during model run.

// Sub-windows knowledge could be removed if vector of obs was used in
// weak constraint 4D-Var. YT

template <typename MODEL, typename STATE> class Observer : public PostBase<STATE> {
  typedef Locations<MODEL>           Locations_;
  typedef ModelAtLocations<MODEL>    GOM_;
  typedef ObsAuxControl<MODEL>       ObsAuxCtrl_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsOperator<MODEL>         ObsOperator_;
  typedef LinearObsOperator<MODEL>   LinearObsOperator_;
  typedef ObservationSpace<MODEL>    ObsSpace_;

 public:
  Observer(const ObsSpace_ &, const ObsOperator_ &,
           const Observations_ &, const ObsAuxCtrl_ &,
           const util::Duration & tslot = util::Duration(0), const bool subwin = false,
           boost::shared_ptr<LinearObsOperator_> htlad = boost::shared_ptr<LinearObsOperator_>() );
  ~Observer() {}

  Observations_ * release() {return yobs_.release();}

 private:
// Methods
  void doInitialize(const STATE &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const STATE &) override;
  void doFinalize(const STATE &) override;

// Obs operator
  ObsSpace_ obspace_;
  const ObsOperator_ hop_;
  boost::shared_ptr<LinearObsOperator_> htlad_;

// Data
  std::auto_ptr<Observations_> yobs_;
  const ObsAuxCtrl_ & ybias_;

  util::DateTime winbgn_;   //!< Begining of assimilation window
  util::DateTime winend_;   //!< End of assimilation window
  util::DateTime bgn_;      //!< Begining of currently active observations
  util::DateTime end_;      //!< End of currently active observations
  util::Duration hslot_;    //!< Half time slot
  const bool subwindows_;

  boost::scoped_ptr<GOM_> gom_;
};

// ====================================================================================

template <typename MODEL, typename STATE>
Observer<MODEL, STATE>::Observer(const ObsSpace_ & obsdb,
                                 const ObsOperator_ & hop,
                                 const Observations_ & yobs,
                                 const ObsAuxCtrl_ & ybias,
                                 const util::Duration & tslot, const bool swin,
                                 boost::shared_ptr<LinearObsOperator_> htlad)
  : PostBase<STATE>(), obspace_(obsdb), hop_(hop), htlad_(htlad),
    yobs_(new Observations_(obsdb)), ybias_(ybias),
    winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()),
    bgn_(winbgn_), end_(winend_), hslot_(tslot/2), subwindows_(swin),
    gom_()
{}
// -----------------------------------------------------------------------------
template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::doInitialize(const STATE & xx,
              const util::DateTime & end, const util::Duration & tstep) {
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
// Pass the Geometry for IFS -- Bad...
  gom_.reset(new GOM_(obspace_, hop_.variables(), bgn_, end_, xx.geometry()));
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::doProcessing(const STATE & xx) {
  util::DateTime t1(xx.validTime()-hslot_);
  util::DateTime t2(xx.validTime()+hslot_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

// Get locations info for interpolator
  Locations_ locs(obspace_, t1, t2);

// Interpolate state variables to obs locations
  xx.interpolate(locs, *gom_);
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename STATE>
void Observer<MODEL, STATE>::doFinalize(const STATE &) {
  if (htlad_) htlad_->setTrajectory(*gom_, ybias_);
  yobs_->runObsOperator(hop_, *gom_, ybias_);
  gom_.reset();
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVER_H_
