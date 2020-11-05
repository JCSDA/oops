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

#include <algorithm>
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
            QCData_ &, const int iteration = 0);
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
  util::Duration hslot_;    //!< Half time slot

  std::vector<std::shared_ptr<Observer_>> observers_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observers<MODEL, OBS>::Observers(const eckit::Configuration & conf, const ObsSpaces_ & obsdb,
                                 const ObsAuxCtrls_ & ybias, QCData_ & qc, const int iteration)
  : PostBase<State_>(),
    obspace_(obsdb), yobs_(obsdb),
    winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()), hslot_(0), observers_(0)
{
  Log::trace() << "Observers::Observers starting" << std::endl;
  std::vector<eckit::LocalConfiguration> typeconf = conf.getSubConfigurations();
  ASSERT(obsdb.size() == typeconf.size());
  observers_.reserve(obsdb.size());
  for (size_t jj = 0; jj < obsdb.size(); ++jj) {
    observers_.emplace_back(new Observer_(typeconf[jj], obsdb[jj],
                       ybias[jj], yobs_[jj], qc.qcFlags(jj), qc.obsErrors(jj), iteration));
  }
  Log::trace() << "Observers::Observers done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::doInitialize(const State_ & xx, const util::DateTime & end,
                                         const util::Duration & tstep) {
  Log::trace() << "Observers::doInitialize start" << std::endl;
  const util::DateTime bgn(xx.validTime());
  hslot_ = tstep/2;

  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    observers_[jj]->doInitialize(xx, winbgn_, winend_);
  }
  Log::trace() << "Observers::doInitialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::doProcessing(const State_ & xx) {
  Log::trace() << "Observers::doProcessing start" << std::endl;
  util::DateTime t1 = std::max(xx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(xx.validTime()+hslot_, winend_);

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
