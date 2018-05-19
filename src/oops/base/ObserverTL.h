/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVERTL_H_
#define OOPS_BASE_OBSERVERTL_H_

#include <memory>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Departures.h"
#include "oops/base/LinearObsOperators.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostBaseTL.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

/// Computes observation equivalent TL and AD to/from increments.

template <typename MODEL, typename INCR> class ObserverTL : public PostBaseTL<INCR> {
  typedef Departures<MODEL>          Departures_;
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef LinearObsOperators<MODEL>  LinearObsOperator_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncr_;
  typedef ObsSpaces<MODEL>           ObsSpace_;

 public:
  ObserverTL(const ObsSpace_ &, const LinearObsOperator_ &, const ObsAuxIncr_ &,
             const util::Duration &, const bool subwin = false);
  ~ObserverTL() {}

  Departures_ * releaseOutputFromTL() override {return ydep_.release();}

 private:
// Methods
  void doInitializeTL(const INCR &, const util::DateTime &, const util::Duration &) override;
  void doProcessingTL(const INCR &) override;
  void doFinalizeTL(const INCR &) override;

// Obs operator
  const ObsSpace_ & obspace_;
  const LinearObsOperator_ & hoptlad_;

// Data
  std::auto_ptr<Departures_> ydep_;
  const ObsAuxIncr_ & ybias_;

  util::DateTime winbgn_;   //!< Begining of assimilation window
  util::DateTime winend_;   //!< End of assimilation window
  util::DateTime bgn_;      //!< Begining of currently active observations
  util::DateTime end_;      //!< End of currently active observations
  util::Duration hslot_;    //!< Half time slot
  const bool subwindows_;

  std::vector<boost::shared_ptr<GeoVaLs_> > gvals_;
};

// ====================================================================================

template <typename MODEL, typename INCR>
ObserverTL<MODEL, INCR>::ObserverTL(const ObsSpace_ & obsdb,
                                    const LinearObsOperator_ & hoptlad,
                                    const ObsAuxIncr_ & ybias,
                                    const util::Duration & tslot, const bool subwin)
  : PostBaseTL<INCR>(obsdb.windowStart(), obsdb.windowEnd()),
    obspace_(obsdb), hoptlad_(hoptlad),
    ydep_(new Departures_(obspace_)), ybias_(ybias),
    winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()),
    bgn_(winbgn_), end_(winend_), hslot_(tslot/2), subwindows_(subwin)
{}
// -----------------------------------------------------------------------------
template <typename MODEL, typename INCR>
void ObserverTL<MODEL, INCR>::doInitializeTL(const INCR & dx,
                 const util::DateTime & end, const util::Duration & tstep) {
  const util::DateTime bgn(dx.validTime());
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

  for (std::size_t jj = 0; jj < obspace_.size(); ++jj) {
    boost::shared_ptr<GeoVaLs_>
      gom(new GeoVaLs_(obspace_[jj].locations(bgn_, end_), hoptlad_.variables(jj)));
    gvals_.push_back(gom);
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename INCR>
void ObserverTL<MODEL, INCR>::doProcessingTL(const INCR & dx) {
  util::DateTime t1(dx.validTime()-hslot_);
  util::DateTime t2(dx.validTime()+hslot_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

// Interpolate state variables to obs locations
  for (std::size_t jj = 0; jj < obspace_.size(); ++jj) {
    dx.interpolateTL(obspace_[jj].locations(t1, t2), hoptlad_.variables(jj), *gvals_.at(jj));
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename INCR>
void ObserverTL<MODEL, INCR>::doFinalizeTL(const INCR &) {
  for (std::size_t jj = 0; jj < obspace_.size(); ++jj) {
    hoptlad_[jj].obsEquivTL(*gvals_.at(jj), (*ydep_)[jj], ybias_);
  }
  gvals_.clear();
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERTL_H_
