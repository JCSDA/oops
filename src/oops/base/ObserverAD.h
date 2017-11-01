/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVERAD_H_
#define OOPS_BASE_OBSERVERAD_H_

#include <memory>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/base/Departures.h"
#include "oops/base/Observations.h"
#include "oops/base/PostBaseAD.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ModelAtLocations.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/LinearObsOperator.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

/// Computes observation equivalent AD and AD to/from increments.

template <typename MODEL, typename INCR> class ObserverAD : public PostBaseAD<INCR> {
  typedef Departures<MODEL>          Departures_;
  typedef Locations<MODEL>           Locations_;
  typedef ModelAtLocations<MODEL>    GOM_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncr_;
  typedef LinearObsOperator<MODEL>   LinearObsOperator_;
  typedef ObservationSpace<MODEL>    ObsSpace_;

 public:
  ObserverAD(const ObsSpace_ &, const LinearObsOperator_ &,
             boost::shared_ptr<const Departures_>, ObsAuxIncr_ &,
             const util::Duration &, const bool subwin = false);
  ~ObserverAD() {}

 private:
// Methods
  void doFirstAD(INCR &, const util::DateTime &, const util::Duration &) override;
  void doProcessingAD(INCR &) override;
  void doLastAD(INCR &) override;

// Obs operator
  ObsSpace_ obspace_;
  const LinearObsOperator_ hoptlad_;

// Data
  boost::shared_ptr<const Departures_> ydep_;
  ObsAuxIncr_ & ybias_;

  util::DateTime winbgn_;   //!< Begining of assimilation window
  util::DateTime winend_;   //!< End of assimilation window
  util::DateTime bgn_;      //!< Begining of currently active observations
  util::DateTime end_;      //!< End of currently active observations
  util::Duration hslot_;    //!< Half time slot
  const bool subwindows_;

  boost::scoped_ptr<GOM_> gom_;
};

// ====================================================================================

template <typename MODEL, typename INCR>
ObserverAD<MODEL, INCR>::ObserverAD(const ObsSpace_ & obsdb,
                                    const LinearObsOperator_ & hoptlad,
                                    boost::shared_ptr<const Departures_> ydep,
                                    ObsAuxIncr_ & ybias,
                                    const util::Duration & tslot, const bool subwin)
  : PostBaseAD<INCR>(obsdb.windowStart(), obsdb.windowEnd()),
    obspace_(obsdb), hoptlad_(hoptlad),
    ydep_(ydep), ybias_(ybias),
    winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()),
    bgn_(winbgn_), end_(winend_), hslot_(tslot/2), subwindows_(subwin)
{}
// -----------------------------------------------------------------------------
template <typename MODEL, typename INCR>
void ObserverAD<MODEL, INCR>::doFirstAD(INCR & dx, const util::DateTime & bgn,
                                        const util::Duration & tstep) {
  if (hslot_ == util::Duration(0)) hslot_ = tstep/2;
  const util::DateTime end(dx.validTime());
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
  gom_.reset(new GOM_(obspace_, hoptlad_.variables(), bgn_, end_, dx.geometry()));
  ydep_->runObsOperatorAD(hoptlad_, *gom_, ybias_);
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename INCR>
void ObserverAD<MODEL, INCR>::doProcessingAD(INCR & dx) {
  util::DateTime t1(dx.validTime()-hslot_);
  util::DateTime t2(dx.validTime()+hslot_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

// Get locations info for interpolator
  Locations_ locs(obspace_, t1, t2);

// Interpolate state variables to obs locations
  dx.interpolateAD(locs, *gom_);
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename INCR>
void ObserverAD<MODEL, INCR>::doLastAD(INCR &) {
  gom_.reset();
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERAD_H_
