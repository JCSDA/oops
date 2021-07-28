/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_GETVALUESPOST_H_
#define OOPS_BASE_GETVALUESPOST_H_

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "oops/base/ObsSpaces.h"
#include "oops/base/PostBase.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/GetValues.h"
#include "oops/interface/Locations.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// \brief Fills GeoVaLs with requested variables at requested locations:
/// - as a postprocessor during model run.
/// - as a method fill() on State4D
template <typename MODEL, typename OBS>
class GetValuesPost : public PostBase<State<MODEL>> {
  typedef GeoVaLs<OBS>              GeoVaLs_;
  typedef Locations<OBS>            Locations_;
  typedef ObsSpaces<OBS>            ObsSpaces_;
  typedef State<MODEL>              State_;
  typedef State4D<MODEL>            State4D_;
  typedef GetValues<MODEL, OBS>     GetValues_;

  typedef std::vector<std::unique_ptr<GetValues_>>  GetValuesVec_;
  typedef std::vector<std::unique_ptr<GeoVaLs_>>    GeoVaLsVec_;
  typedef std::vector<std::unique_ptr<Locations_>>  LocationsVec_;

 public:
/// \brief Saves Locations and Variables to be processed
  GetValuesPost(const ObsSpaces_ &,
                const LocationsVec_ &, const std::vector<Variables> &,
                const std::vector<eckit::LocalConfiguration> &);

/// \brief Returns geovals filled in during the model run
  const GeoVaLsVec_ & geovals() const {return geovals_;}

/// \brief fills in GeoVaLs looping through State4D
  void fill(const State4D_ &);

 private:
/// \brief initialization before model run: sets up GetValues and allocate GeoVaLs
  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
/// \brief called at each model step: fill in GeoVaLs for the current time slot
  void doProcessing(const State_ &) override;

// Data
  util::DateTime winbgn_;   /// Begining of assimilation window
  util::DateTime winend_;   /// End of assimilation window
  util::Duration hslot_;    /// Half time slot

  const LocationsVec_ & locations_;   /// locations of observations
  const std::vector<Variables> geovars_;       /// Variables needed from model
  GetValuesVec_ getvals_;             /// GetValues used to fill in GeoVaLs
  GeoVaLsVec_ geovals_;               /// GeoVaLs that are filled in
  const std::vector<eckit::LocalConfiguration> getvalsconfs_;   /// configuration object
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GetValuesPost<MODEL, OBS>::GetValuesPost(const ObsSpaces_ & obsdb,
                                         const LocationsVec_ & locations,
                                         const std::vector<Variables> & vars,
                                         const std::vector<eckit::LocalConfiguration> & confs)
  : PostBase<State_>(),
    winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()), hslot_(),
    locations_(locations), geovars_(vars), getvalsconfs_(confs)
{
  Log::trace() << "GetValuesPost::GetValuesPost" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void GetValuesPost<MODEL, OBS>::fill(const State4D_ & xx) {
  Log::trace() << "GetValuesPost::fill start" << std::endl;
  const size_t nstates = xx.size();
  util::Duration tstep = winend_ - winbgn_;   // for a single state
  // if using several states, compute the timestep and check that it's the same
  // for all states
  if (nstates > 1) {
    tstep = xx[1].validTime() - xx[0].validTime();
    for (size_t ii = 1; ii < nstates; ++ii) {
      ASSERT(tstep == (xx[ii].validTime() - xx[ii-1].validTime()));
    }
  }
  // run GetValues postprocessor looping through all the states
  doInitialize(xx[0], xx[nstates-1].validTime(), tstep);
  for (size_t ii = 0; ii < nstates; ++ii) {
    doProcessing(xx[ii]);
  }
  Log::trace() << "GetValuesPost::fill done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuesPost<MODEL, OBS>::doInitialize(const State_ & xx, const util::DateTime & end,
                                             const util::Duration & tstep) {
  Log::trace() << "GetValuesPost::doInitialize start" << std::endl;
  hslot_ = tstep/2;

  for (size_t jj = 0; jj < locations_.size(); ++jj) {
    getvals_.emplace_back(new GetValues_(xx.geometry(), *locations_[jj], getvalsconfs_[jj]));
    geovals_.emplace_back(new GeoVaLs_(*locations_[jj], geovars_[jj],
                                       xx.geometry().variableSizes(geovars_[jj])));
  }

  Log::trace() << "GetValuesPost::doInitialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuesPost<MODEL, OBS>::doProcessing(const State_ & xx) {
  Log::trace() << "GetValuesPost::doProcessing start" << std::endl;
  util::DateTime t1 = std::max(xx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(xx.validTime()+hslot_, winend_);

// Get state variables at obs locations
  for (size_t jj = 0; jj < getvals_.size(); ++jj) {
    getvals_[jj]->fillGeoVaLs(xx, t1, t2, *geovals_[jj]);
  }
  Log::trace() << "GetValuesPost::doProcessing done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_GETVALUESPOST_H_
