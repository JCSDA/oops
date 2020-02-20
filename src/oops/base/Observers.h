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
#include "oops/base/Variables.h"
#include "oops/interface/ObsDataVector.h"
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
class Observers : public PostBase<STATE>,
                  public util::Printable {
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObsAuxControls<MODEL>      ObsAuxCtrls_;
  typedef Observations<MODEL>        Observations_;
  typedef Observer<MODEL>            Observer_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  template <typename DATA> using ObsDataPtr_ = boost::shared_ptr<ObsDataVector<MODEL, DATA> >;
  template <typename DATA> using ObsDataVectors_ = std::vector<ObsDataPtr_<DATA> >;
 public:
  Observers(const eckit::Configuration &, const ObsSpaces_ & obsdb, const ObsAuxCtrls_ &,
            ObsDataVectors_<int> = ObsDataVectors_<int>(),
            ObsDataVectors_<float> = ObsDataVectors_<float>(),
            const util::Duration & tslot = util::Duration(0), const bool subwin = false);
  ~Observers() {}

  const Observations_ & hofx() {return yobs_;}

 private:
// Methods
  void doInitialize(const STATE &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const STATE &) override;
  void doFinalize(const STATE &) override;
  void print(std::ostream &) const override;

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

template <typename MODEL, typename STATE>
Observers<MODEL, STATE>::Observers(const eckit::Configuration & conf,
                                   const ObsSpaces_ & obsdb,
                                   const ObsAuxCtrls_ & ybias,
                                   ObsDataVectors_<int> qcflags,
                                   ObsDataVectors_<float> obserr,
                                   const util::Duration & tslot, const bool swin)
  : PostBase<STATE>(),
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
  if (qcflags.size() == obsdb.size() && obserr.size() == obsdb.size()) {
    for (size_t jj = 0; jj < obsdb.size(); ++jj) {
      typeconf[jj].set("iteration", iterout);
      observers_.emplace_back(new Observer_(typeconf[jj], obsdb[jj],
                                  ybias[jj], yobs_[jj], qcflags[jj], obserr[jj]));
    }
  } else if (qcflags.size() == 0 && obserr.size() == 0) {
    for (size_t jj = 0; jj < obsdb.size(); ++jj) {
      typeconf[jj].set("iteration", iterout);
      observers_.emplace_back(new Observer_(typeconf[jj], obsdb[jj],
                                  ybias[jj], yobs_[jj]));
    }
  } else {
    ABORT("Observers: have to provide qcflags and obserrs for all or none of the ObsTypes");
  }

  Log::trace() << "Observers::Observers done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
void Observers<MODEL, STATE>::doInitialize(const STATE & xx,
                                          const util::DateTime & end,
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

template <typename MODEL, typename STATE>
void Observers<MODEL, STATE>::doProcessing(const STATE & xx) {
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

template <typename MODEL, typename STATE>
void Observers<MODEL, STATE>::doFinalize(const STATE &) {
  Log::trace() << "Observers::doFinalize start" << std::endl;
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    observers_[jj]->doFinalize();
  }
  Log::trace() << "Observers::doFinalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename STATE>
void Observers<MODEL, STATE>::print(std::ostream &) const {}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERS_H_
