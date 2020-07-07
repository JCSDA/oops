/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVERS_H_
#define OOPS_BASE_OBSERVERS_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/ObsAuxControls.h"
#include "oops/base/Observations.h"
#include "oops/base/Observer.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostBase.h"
#include "oops/base/QCData.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// Computes observation equivalent during model run.

// Sub-windows knowledge could be removed if vector of obs was used in
// weak constraint 4D-Var. YT

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
class Observers : public PostBase<State<MODEL>> {
  typedef GeoVaLs<OBS>             GeoVaLs_;
  typedef ObsAuxControls<OBS>      ObsAuxCtrls_;
  typedef Observations<OBS>        Observations_;
  typedef Observer<MODEL, OBS>     Observer_;
  typedef ObsSpaces<OBS>           ObsSpaces_;
  typedef QCData<OBS>              QCData_;
  typedef State<MODEL>               State_;

 public:
  Observers(const eckit::Configuration &, const ObsSpaces_ & obsdb, const ObsAuxCtrls_ &,
            QCData_ &,
            const util::Duration & tslot = util::Duration(0), const bool subwin = false);
  ~Observers() {}

  const Observations_ & hofx() {return yobs_;}

 private:
// Methods
  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const State_ &) override;
  void doFinalize(const State_ &) override;

// Data
  ObsSpaces_ obspace_;
  Observations_ yobs_;

  util::DateTime winbgn_;   //!< Begining of assimilation window
  util::DateTime winend_;   //!< End of assimilation window
  util::DateTime bgn_;      //!< Begining of currently active observations
  util::DateTime end_;      //!< End of currently active observations
  util::Duration hslot_;    //!< Half time slot
  const bool subwindows_;

  std::vector<std::shared_ptr<Observer_>> observers_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observers<MODEL, OBS>::Observers(const eckit::Configuration & conf, const ObsSpaces_ & obsdb,
                                 const ObsAuxCtrls_ & ybias, QCData_ & qc,
                                 const util::Duration & tslot, const bool swin)
  : PostBase<State_>(),
    obspace_(obsdb), yobs_(obsdb),
    winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()),
    bgn_(winbgn_), end_(winend_), hslot_(tslot/2), subwindows_(swin),
    observers_(0)
{
  Log::trace() << "Observers::Observers starting" << std::endl;

  const int iterout = conf.getInt("iteration", 0);
  std::vector<eckit::LocalConfiguration> typeconf;
  conf.get("ObsTypes", typeconf);
  observers_.reserve(obsdb.size());
  for (size_t jj = 0; jj < obsdb.size(); ++jj) {
    typeconf[jj].set("iteration", iterout);
    observers_.emplace_back(new Observer_(typeconf[jj], obsdb[jj],
                                ybias[jj], yobs_[jj], qc.qcFlags(jj), qc.obsErrors(jj)));
  }
  Log::trace() << "Observers::Observers done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::doInitialize(const State_ & xx, const util::DateTime & end,
                                         const util::Duration & tstep) {
  Log::trace() << "Observers::doInitialize start" << std::endl;
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

  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    observers_[jj]->doInitialize(xx, bgn_, end_);
  }
  Log::trace() << "Observers::doInitialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::doProcessing(const State_ & xx) {
  Log::trace() << "Observers::doProcessing start" << std::endl;
  util::DateTime t1(xx.validTime()-hslot_);
  util::DateTime t2(xx.validTime()+hslot_);
  if (t1 < bgn_) t1 = bgn_;
  if (t2 > end_) t2 = end_;

// Get state variables at obs locations
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    observers_[jj]->doProcessing(xx, t1, t2);
  }
  Log::trace() << "Observers::doProcessing done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::doFinalize(const State_ &) {
  Log::trace() << "Observers::doFinalize start" << std::endl;
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    observers_[jj]->doFinalize();
  }
  Log::trace() << "Observers::doFinalize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERS_H_
