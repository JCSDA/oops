/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVERTLAD_H_
#define OOPS_BASE_OBSERVERTLAD_H_

#include <memory>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/base/Departures.h"
#include "oops/base/InterpolatorTLAD.h"
#include "oops/base/LinearObsOperators.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsFilters.h"
#include "oops/base/ObsOperators.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Increment.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

/// Computes observation equivalent TL and AD to/from increments.

template <typename MODEL>
class ObserverTLAD : public PostBaseTLAD<MODEL> {
  typedef Departures<MODEL>          Departures_;
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef Increment<MODEL>           Increment_;
  typedef InterpolatorTLAD<MODEL>    InterpolatorTLAD_;
  typedef LinearObsOperators<MODEL>  LinearObsOperators_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsAuxControl<MODEL>       ObsAuxCtrl_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncr_;
  typedef ObsFilters<MODEL>          ObsFilters_;
  typedef ObsOperators<MODEL>        ObsOperators_;
  typedef ObsSpaces<MODEL>           ObsSpace_;
  typedef State<MODEL>               State_;

 public:
  ObserverTLAD(const ObsSpace_ &, const ObsOperators_ &, const ObsAuxCtrl_ &, const ObsFilters_ &,
               const util::Duration & tslot = util::Duration(0), const bool subwin = false);
  ~ObserverTLAD() {}

  Observations_ * release() {return observer_.release();}
  Departures_ * releaseOutputFromTL() override {return ydeptl_.release();}
  void setupTL(const ObsAuxIncr_ &);
  void setupAD(boost::shared_ptr<const Departures_>, ObsAuxIncr_ &);

 private:
// Methods
  void doInitializeTraj(const State_ &,
                        const util::DateTime &, const util::Duration &) override;
  void doProcessingTraj(const State_ &) override;
  void doFinalizeTraj(const State_ &) override;

  void doInitializeTL(const Increment_ &,
                      const util::DateTime &, const util::Duration &) override;
  void doProcessingTL(const Increment_ &) override;
  void doFinalizeTL(const Increment_ &) override;

  void doFirstAD(Increment_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingAD(Increment_ &) override;
  void doLastAD(Increment_ &) override;

// Obs operator
  const ObsSpace_ & obspace_;
  LinearObsOperators_ hoptlad_;
  Observer<MODEL, State_> observer_;

// Data
  std::auto_ptr<Departures_> ydeptl_;
  const ObsAuxIncr_ * ybiastl_;
  boost::shared_ptr<const Departures_> ydepad_;
  ObsAuxIncr_ * ybiasad_;
  InterpolatorTLAD_ interptlad_;

  util::DateTime winbgn_;   //!< Begining of assimilation window
  util::DateTime winend_;   //!< End of assimilation window
  util::DateTime bgn_;      //!< Begining of currently active observations
  util::DateTime end_;      //!< End of currently active observations
  util::Duration hslot_;    //!< Half time slot
  const bool subwindows_;

  std::vector<boost::shared_ptr<GeoVaLs_> > gvals_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObserverTLAD<MODEL>::ObserverTLAD(const ObsSpace_ & obsdb,
                                  const ObsOperators_ & hop,
                                  const ObsAuxCtrl_ & ybias,
                                  const ObsFilters_ & filters,
                                  const util::Duration & tslot, const bool subwin)
  : PostBaseTLAD<MODEL>(obsdb.windowStart(), obsdb.windowEnd()),
    obspace_(obsdb), hoptlad_(obspace_),
    observer_(obspace_, hop, ybias, filters, tslot, subwin),
    ydeptl_(), ybiastl_(), ydepad_(), ybiasad_(), interptlad_(obspace_.size()),
    winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()),
    bgn_(winbgn_), end_(winend_), hslot_(tslot/2), subwindows_(subwin)
{
  Log::trace() << "ObserverTLAD::ObserverTLAD" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doInitializeTraj(const State_ & xx,
               const util::DateTime & end, const util::Duration & tstep) {
  Log::trace() << "ObserverTLAD::doInitializeTraj start" << std::endl;
  observer_.initialize(xx, end, tstep);
  Log::trace() << "ObserverTLAD::doInitializeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doProcessingTraj(const State_ & xx) {
  Log::trace() << "ObserverTLAD::doProcessingTraj start" << std::endl;
  observer_.processTraj(xx, interptlad_);
  Log::trace() << "ObserverTLAD::doProcessingTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doFinalizeTraj(const State_ & xx) {
  Log::trace() << "ObserverTLAD::doFinalizeTraj start" << std::endl;
  observer_.finalizeTraj(xx, hoptlad_);
  Log::trace() << "ObserverTLAD::doFinalizeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::setupTL(const ObsAuxIncr_ & ybias) {
  Log::trace() << "ObserverTLAD::setupTL start" << std::endl;
  ydeptl_.reset(new Departures_(obspace_));
  ybiastl_ = &ybias;
  Log::trace() << "ObserverTLAD::setupTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doInitializeTL(const Increment_ & dx,
                   const util::DateTime & end, const util::Duration & tstep) {
  Log::trace() << "ObserverTLAD::doInitializeTL start" << std::endl;
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
  Log::trace() << "ObserverTLAD::doInitializeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doProcessingTL(const Increment_ & dx) {
  Log::trace() << "ObserverTLAD::doProcessingTL start" << std::endl;
  util::DateTime t1(dx.validTime()-hslot_);
  util::DateTime t2(dx.validTime()+hslot_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

// Get increment variables at obs locations
  for (std::size_t jj = 0; jj < obspace_.size(); ++jj) {
    dx.getValuesTL(obspace_[jj].locations(t1, t2), hoptlad_.variables(jj),
                   *gvals_.at(jj), interptlad_[jj]);
  }
  Log::trace() << "ObserverTLAD::doProcessingTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doFinalizeTL(const Increment_ &) {
  Log::trace() << "ObserverTLAD::doFinalizeTL start" << std::endl;
  for (std::size_t jj = 0; jj < obspace_.size(); ++jj) {
    hoptlad_[jj].obsEquivTL(*gvals_.at(jj), (*ydeptl_)[jj], *ybiastl_);
  }
  gvals_.clear();
  Log::trace() << "ObserverTLAD::doFinalizeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::setupAD(boost::shared_ptr<const Departures_> ydep,
                                  ObsAuxIncr_ & ybias) {
  Log::trace() << "ObserverTLAD::setupAD start" << std::endl;
  ydepad_ = ydep;
  ybiasad_ = &ybias;
  Log::trace() << "ObserverTLAD::setupAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doFirstAD(Increment_ & dx, const util::DateTime & bgn,
                                    const util::Duration & tstep) {
  Log::trace() << "ObserverTLAD::doFirstAD start" << std::endl;
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

  for (std::size_t jj = 0; jj < obspace_.size(); ++jj) {
    boost::shared_ptr<GeoVaLs_>
      gom(new GeoVaLs_(obspace_[jj].locations(bgn_, end_), hoptlad_.variables(jj)));
    hoptlad_[jj].obsEquivAD(*gom, (*ydepad_)[jj], *ybiasad_);
    gvals_.push_back(gom);
  }
  Log::trace() << "ObserverTLAD::doFirstAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doProcessingAD(Increment_ & dx) {
  Log::trace() << "ObserverTLAD::doProcessingAD start" << std::endl;
  util::DateTime t1(dx.validTime()-hslot_);
  util::DateTime t2(dx.validTime()+hslot_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

// Adjoint of get increment variables at obs locations
  for (std::size_t jj = 0; jj < obspace_.size(); ++jj) {
    dx.getValuesAD(obspace_[jj].locations(t1, t2), hoptlad_.variables(jj),
                   *gvals_.at(jj), interptlad_[jj]);
  }
  Log::trace() << "ObserverTLAD::doProcessingAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doLastAD(Increment_ &) {
  Log::trace() << "ObserverTLAD::doLastAD start" << std::endl;
  gvals_.clear();
  Log::trace() << "ObserverTLAD::doLastAD done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERTLAD_H_
