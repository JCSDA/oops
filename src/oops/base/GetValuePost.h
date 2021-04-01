/*
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_GETVALUEPOST_H_
#define OOPS_BASE_GETVALUEPOST_H_

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/PostBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/GetValues.h"
#include "oops/interface/Locations.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// \brief Fills GeoVaLs with requested variables at requested locations during model run
template <typename MODEL, typename OBS>
class GetValuePost : public PostBase<State<MODEL>> {
  typedef Geometry<MODEL>           Geometry_;
  typedef GeoVaLs<OBS>              GeoVaLs_;
  typedef GetValues<MODEL, OBS>     GetValues_;
  typedef Locations<OBS>            Locations_;
  typedef State<MODEL>              State_;

 public:
/// \brief Saves Locations and Variables to be processed
  GetValuePost(const eckit::Configuration &, const Geometry_ &,
               const util::DateTime &, const util::DateTime &,
               const Locations_ &, const Variables &);

/// \brief Returns geovals filled in during the model run
  std::unique_ptr<GeoVaLs_> releaseGeoVaLs();

 private:
/// \brief initialization before model run: sets up GetValues and allocate GeoVaLs
  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
/// \brief called at each model step: fill in GeoVaLs for the current time slot
  void doProcessing(const State_ &) override;

// Data
  util::DateTime winbgn_;   /// Begining of assimilation window
  util::DateTime winend_;   /// End of assimilation window
  util::Duration hslot_;    /// Half time slot

  const Locations_ & locations_;       /// locations of observations
  const Variables geovars_;            /// Variables needed from model
  GetValues_ getvals_;                 /// GetValues used to fill in GeoVaLs
  std::unique_ptr<GeoVaLs_> geovals_;  /// GeoVaLs that are filled in
  bool initialized_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GetValuePost<MODEL, OBS>::GetValuePost(const eckit::Configuration & conf, const Geometry_ & geom,
                                       const util::DateTime & bgn, const util::DateTime & end,
                                       const Locations_ & locations, const Variables & vars)
  : PostBase<State_>(),
    winbgn_(bgn), winend_(end), hslot_(), locations_(locations), geovars_(vars),
    getvals_(geom, locations_, conf), geovals_(), initialized_(false)
{
  Log::trace() << "GetValuePost::GetValuePost" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuePost<MODEL, OBS>::doInitialize(const State_ & xx, const util::DateTime &,
                                            const util::Duration & tstep) {
  Log::trace() << "GetValuePost::doInitialize start" << std::endl;
  hslot_ = tstep/2;

  geovals_.reset(new GeoVaLs_(locations_, geovars_));

  initialized_ = true;
  Log::trace() << "GetValuePost::doInitialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuePost<MODEL, OBS>::doProcessing(const State_ & xx) {
  Log::trace() << "GetValuePost::doProcessing start" << std::endl;
  ASSERT(initialized_);
  util::DateTime t1 = std::max(xx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(xx.validTime()+hslot_, winend_);

// Get state variables at obs locations
  getvals_.fillGeoVaLs(xx, t1, t2, *geovals_);
  Log::trace() << "GetValuePost::doProcessing done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<GeoVaLs<OBS>> GetValuePost<MODEL, OBS>::releaseGeoVaLs() {
  Log::trace() << "GetValuePost::releaseGeoVaLs" << std::endl;
  initialized_ = false;
  // Release ownership of GeoVaLs
  return std::move(geovals_);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_GETVALUEPOST_H_
