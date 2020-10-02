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

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Departures.h"
#include "oops/base/GeneralizedDepartures.h"
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

template <typename MODEL, typename OBS>
class ObserversTLAD : public PostBaseTLAD<MODEL> {
  typedef Departures<OBS>           Departures_;
  typedef GeoVaLs<OBS>              GeoVaLs_;
  typedef Increment<MODEL>          Increment_;
  typedef Observations<OBS>         Observations_;
  typedef ObsAuxControls<OBS>       ObsAuxCtrls_;
  typedef ObsAuxIncrements<OBS>     ObsAuxIncrs_;
  typedef ObserverTLAD<MODEL, OBS>  ObserverTLAD_;
  typedef ObsSpaces<OBS>            ObsSpaces_;
  typedef State<MODEL>              State_;

 public:
  ObserversTLAD(const eckit::Configuration &,
                const ObsSpaces_ &, const ObsAuxCtrls_ &,
                const util::Duration & tslot = util::Duration(0), const bool subwin = false);
  ~ObserversTLAD() {}

  std::unique_ptr<GeneralizedDepartures> releaseOutputFromTL() override {return std::move(ydeptl_);}
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
  std::unique_ptr<Departures_> ydeptl_;
  const ObsAuxIncrs_ * ybiastl_;
  std::shared_ptr<const Departures_> ydepad_;
  ObsAuxIncrs_ * ybiasad_;

  util::DateTime winbgn_;   //!< Begining of assimilation window
  util::DateTime winend_;   //!< End of assimilation window
  util::Duration hslot_, hslottraj_;    //!< Half time slot
};

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
ObserversTLAD<MODEL, OBS>::ObserversTLAD(const eckit::Configuration & config,
                                  const ObsSpaces_ & obsdb,
                                  const ObsAuxCtrls_ & ybias,
                                  const util::Duration & tslot, const bool subwin)
  : PostBaseTLAD<MODEL>(obsdb.windowStart(), obsdb.windowEnd()),
    observerstlad_(), obspace_(obsdb),
    ydeptl_(), ybiastl_(), ydepad_(), ybiasad_(),
    winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()),
    hslot_(tslot/2), hslottraj_(tslot/2)
{
  // setup observers
  std::vector<eckit::LocalConfiguration> typeconf;
  config.get("observations", typeconf);
  for (std::size_t jobs = 0; jobs < obsdb.size(); ++jobs) {
    // Set LinearObsOperator section to ObsOperator section if not available
    if (!typeconf[jobs].has("linear obs operator")) {
      typeconf[jobs].set("linear obs operator", typeconf[jobs].getSubConfiguration("obs operator"));
    }
    std::shared_ptr<ObserverTLAD_> tmp(new ObserverTLAD_(typeconf[jobs], obsdb[jobs], ybias[jobs]));
    observerstlad_.push_back(tmp);
  }
  Log::trace() << "ObserversTLAD::ObserversTLAD" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::doInitializeTraj(const State_ & xx,
                                                 const util::DateTime & end,
                                                 const util::Duration & tstep) {
  Log::trace() << "ObserversTLAD::doInitializeTraj start" << std::endl;
// Create full trajectory object

  if (hslottraj_ == util::Duration(0)) hslottraj_ = tstep/2;

  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doInitializeTraj(xx, winbgn_, winend_);
  }
  Log::trace() << "ObserversTLAD::doInitializeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::doProcessingTraj(const State_ & xx) {
  Log::trace() << "ObserversTLAD::doProcessingTraj start" << std::endl;
  util::DateTime t1 = std::max(xx.validTime()-hslottraj_, winbgn_);
  util::DateTime t2 = std::min(xx.validTime()+hslottraj_, winend_);

  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doProcessingTraj(xx, t1, t2);
  }
  Log::trace() << "ObserversTLAD::doProcessingTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::doFinalizeTraj(const State_ & xx) {
  Log::trace() << "ObserversTLAD::doFinalizeTraj start" << std::endl;
  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doFinalizeTraj(xx);
  }
  Log::trace() << "ObserversTLAD::doFinalizeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::setupTL(const ObsAuxIncrs_ & ybiastl) {
  Log::trace() << "ObserversTLAD::setupTL start" << std::endl;
  ydeptl_.reset(new Departures_(obspace_));
  ybiastl_ = &ybiastl;
  Log::trace() << "ObserversTLAD::setupTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::doInitializeTL(const Increment_ & dx,
                                               const util::DateTime & end,
                                               const util::Duration & tstep) {
  Log::trace() << "ObserversTLAD::doInitializeTL start" << std::endl;
  if (hslot_ == util::Duration(0)) hslot_ = tstep/2;
  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doInitializeTL(dx, winbgn_, winend_);
  }
  Log::trace() << "ObserversTLAD::doInitializeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::doProcessingTL(const Increment_ & dx) {
  Log::trace() << "ObserversTLAD::doProcessingTL start" << std::endl;
  util::DateTime t1 = std::max(dx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(dx.validTime()+hslot_, winend_);

  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doProcessingTL(dx, t1, t2);
  }
  Log::trace() << "ObserversTLAD::doProcessingTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::doFinalizeTL(const Increment_ & dx) {
  Log::trace() << "ObserversTLAD::doFinalizeTL start" << std::endl;
  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doFinalizeTL(dx, (*ydeptl_)[jj], (*ybiastl_)[jj]);
  }
  Log::trace() << "ObserversTLAD::doFinalizeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::setupAD(std::shared_ptr<const Departures_> ydepad,
                                        ObsAuxIncrs_ & ybiasad) {
  Log::trace() << "ObserversTLAD::setupAD start" << std::endl;
  ydepad_  = ydepad;
  ybiasad_ = &ybiasad;
  Log::trace() << "ObserversTLAD::setupAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::doFirstAD(Increment_ & dx, const util::DateTime & bgn,
                                          const util::Duration & tstep) {
  Log::trace() << "ObserversTLAD::doFirstAD start" << std::endl;
  if (hslot_ == util::Duration(0)) hslot_ = tstep/2;

  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doFirstAD(dx, (*ydepad_)[jj], (*ybiasad_)[jj], winbgn_, winend_);
  }
  Log::trace() << "ObserversTLAD::doFirstAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::doProcessingAD(Increment_ & dx) {
  Log::trace() << "ObserversTLAD::doProcessingAD start" << std::endl;
  util::DateTime t1 = std::max(dx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(dx.validTime()+hslot_, winend_);

  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doProcessingAD(dx, t1, t2);
  }
  Log::trace() << "ObserversTLAD::doProcessingAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::doLastAD(Increment_ & dx) {
  Log::trace() << "ObserversTLAD::doLastAD start" << std::endl;
  for (std::size_t jj = 0; jj < observerstlad_.size(); ++jj) {
    observerstlad_[jj]->doLastAD(dx);
  }
  Log::trace() << "ObserversTLAD::doLastAD done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERSTLAD_H_
