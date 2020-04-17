/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVERSTLAD_H_
#define OOPS_BASE_OBSERVERSTLAD_H_

#include <memory>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Departures.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsAuxIncrements.h"
#include "oops/base/Observations.h"
#include "oops/base/ObserverTLAD.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

/// Computes observation equivalent TL and AD to/from increments.

template <typename MODEL>
class ObserversTLAD : public PostBaseTLAD<MODEL> {
  typedef Departures<MODEL>          Departures_;
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef Increment<MODEL>           Increment_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsAuxControls<MODEL>      ObsAuxCtrls_;
  typedef ObsAuxIncrements<MODEL>    ObsAuxIncrs_;
  typedef ObserverTLAD<MODEL>        ObserverTLAD_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef State<MODEL>               State_;

 public:
  ObserversTLAD(const eckit::Configuration &,
                const ObsSpaces_ &, const ObsAuxCtrls_ &,
                const util::Duration & tslot = util::Duration(0), const bool subwin = false);
  ~ObserversTLAD() {}

  Observations_ * release() {return yobs_.release();}
  Departures_ * releaseOutputFromTL() override {return ydeptl_.release();}
  void setupTL(const ObsAuxIncrs_ &);
  void setupAD(std::shared_ptr<const Departures_>, ObsAuxIncrs_ &);

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
  std::vector<std::shared_ptr<ObserverTLAD_>> observerstlad_;

// Data
  ObsSpaces_ obspace_;
  std::unique_ptr<Observations_> yobs_;
  std::unique_ptr<Departures_> ydeptl_;
  const ObsAuxIncrs_ * ybiastl_;
  std::shared_ptr<const Departures_> ydepad_;
  ObsAuxIncrs_ * ybiasad_;

  util::DateTime winbgn_;   //!< Begining of assimilation window
  util::DateTime winend_;   //!< End of assimilation window
  util::DateTime bgn_;      //!< Begining of currently active observations
  util::DateTime end_;      //!< End of currently active observations
  util::Duration hslot_, hslottraj_;    //!< Half time slot
  const bool subwindows_;
  unsigned int nwindows_;   //!< number of subwindows (default 1)
  util::Duration winlen_;   //!< length of subwindow (default winend_-winbgn_)
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObserversTLAD<MODEL>::ObserversTLAD(const eckit::Configuration & config,
                                  const ObsSpaces_ & obsdb,
                                  const ObsAuxCtrls_ & ybias,
                                  const util::Duration & tslot, const bool subwin)
  : PostBaseTLAD<MODEL>(obsdb.windowStart(), obsdb.windowEnd()),
    observerstlad_(), obspace_(obsdb),
    yobs_(new Observations_(obspace_)),
    ydeptl_(), ybiastl_(), ydepad_(), ybiasad_(),
    winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()),
    bgn_(winbgn_), end_(winend_), hslot_(tslot/2), hslottraj_(tslot/2),
    subwindows_(subwin), nwindows_(1), winlen_(winend_-winbgn_)
{
  // setup observers
  std::vector<eckit::LocalConfiguration> typeconf;
  config.get("ObsTypes", typeconf);
  for (std::size_t jobs = 0; jobs < obsdb.size(); ++jobs) {
    // Set LinearObsOperator section to ObsOperator section if not available
    if (!typeconf[jobs].has("LinearObsOperator")) {
      typeconf[jobs].set("LinearObsOperator", typeconf[jobs].getSubConfiguration("ObsOperator"));
    }
    std::shared_ptr<ObserverTLAD_> tmp(new ObserverTLAD_(typeconf[jobs], obsdb[jobs],
                                       ybias[jobs], (*yobs_)[jobs]));
    observerstlad_.push_back(tmp);
  }
  Log::trace() << "ObserversTLAD::ObserversTLAD" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserversTLAD<MODEL>::doInitializeTraj(const State_ & xx,
               const util::DateTime & end, const util::Duration & tstep) {
  Log::trace() << "ObserversTLAD::doInitializeTraj start" << std::endl;
// Create full trajectory object

  const util::DateTime bgn(xx.validTime());
  if (hslottraj_ == util::Duration(0)) hslottraj_ = tstep/2;
  if (subwindows_) {
    if (bgn == end) {
      bgn_ = bgn - hslottraj_;
      end_ = end + hslottraj_;
    } else {
      bgn_ = bgn;
      end_ = end;
    }
    winlen_ = end_ - bgn_;
    nwindows_ = (winend_-winbgn_).toSeconds() / winlen_.toSeconds();
  }
  if (bgn_ < winbgn_) bgn_ = winbgn_;
  if (end_ > winend_) end_ = winend_;

  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doInitializeTraj(xx, bgn_, winlen_, nwindows_);
  }
  Log::trace() << "ObserversTLAD::doInitializeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserversTLAD<MODEL>::doProcessingTraj(const State_ & xx) {
  Log::trace() << "ObserversTLAD::doProcessingTraj start" << std::endl;
  util::DateTime t1(xx.validTime()-hslottraj_);
  util::DateTime t2(xx.validTime()+hslottraj_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

  int iwin = (t1-winbgn_).toSeconds() / winlen_.toSeconds();

  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doProcessingTraj(xx, t1, t2, iwin);
  }
  Log::trace() << "ObserversTLAD::doProcessingTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserversTLAD<MODEL>::doFinalizeTraj(const State_ & xx) {
  Log::trace() << "ObserversTLAD::doFinalizeTraj start" << std::endl;
  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doFinalizeTraj(xx);
  }
  Log::trace() << "ObserversTLAD::doFinalizeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserversTLAD<MODEL>::setupTL(const ObsAuxIncrs_ & ybiastl) {
  Log::trace() << "ObserversTLAD::setupTL start" << std::endl;
  ydeptl_.reset(new Departures_(obspace_));
  ybiastl_ = &ybiastl;
  Log::trace() << "ObserversTLAD::setupTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserversTLAD<MODEL>::doInitializeTL(const Increment_ & dx,
                   const util::DateTime & end, const util::Duration & tstep) {
  Log::trace() << "ObserversTLAD::doInitializeTL start" << std::endl;
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

  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doInitializeTL(dx, bgn_, end_);
  }
  Log::trace() << "ObserversTLAD::doInitializeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserversTLAD<MODEL>::doProcessingTL(const Increment_ & dx) {
  Log::trace() << "ObserversTLAD::doProcessingTL start" << std::endl;
  util::DateTime t1(dx.validTime()-hslot_);
  util::DateTime t2(dx.validTime()+hslot_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

// Index for current bin
  unsigned int iwin = (t1 - winbgn_).toSeconds() / winlen_.toSeconds();
  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doProcessingTL(dx, t1, t2, iwin);
  }
  Log::trace() << "ObserversTLAD::doProcessingTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserversTLAD<MODEL>::doFinalizeTL(const Increment_ & dx) {
  Log::trace() << "ObserversTLAD::doFinalizeTL start" << std::endl;
  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doFinalizeTL(dx, (*ydeptl_)[jj], (*ybiastl_)[jj]);
  }
  Log::trace() << "ObserversTLAD::doFinalizeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserversTLAD<MODEL>::setupAD(std::shared_ptr<const Departures_> ydepad,
                                   ObsAuxIncrs_ & ybiasad) {
  Log::trace() << "ObserversTLAD::setupAD start" << std::endl;
  ydepad_  = ydepad;
  ybiasad_ = &ybiasad;
  Log::trace() << "ObserversTLAD::setupAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserversTLAD<MODEL>::doFirstAD(Increment_ & dx, const util::DateTime & bgn,
                                    const util::Duration & tstep) {
  Log::trace() << "ObserversTLAD::doFirstAD start" << std::endl;
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

  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doFirstAD(dx, (*ydepad_)[jj], (*ybiasad_)[jj], bgn_, end_);
  }
  Log::trace() << "ObserversTLAD::doFirstAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserversTLAD<MODEL>::doProcessingAD(Increment_ & dx) {
  Log::trace() << "ObserversTLAD::doProcessingAD start" << std::endl;
  util::DateTime t1(dx.validTime()-hslot_);
  util::DateTime t2(dx.validTime()+hslot_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

// Index for current window
  unsigned int iwin = (t1 - winbgn_).toSeconds() / winlen_.toSeconds();
  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doProcessingAD(dx, t1, t2, iwin);
  }
  Log::trace() << "ObserversTLAD::doProcessingAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserversTLAD<MODEL>::doLastAD(Increment_ & dx) {
  Log::trace() << "ObserversTLAD::doLastAD start" << std::endl;
  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doLastAD(dx);
  }
  Log::trace() << "ObserversTLAD::doLastAD done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERSTLAD_H_
