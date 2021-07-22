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

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"
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
class GetValuePost {
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

/// \brief initialization before model run: sets up GetValues and allocate GeoVaLs
  void initialize(const util::Duration &);
/// \brief called at each model step: fill in GeoVaLs for the current time slot
  void process(const State_ &);

/// Variables that will be required from the State
  const Variables & requiredVariables() const {return geovars_;}

 private:
  util::DateTime winbgn_;   /// Begining of assimilation window
  util::DateTime winend_;   /// End of assimilation window
  util::Duration hslot_;    /// Half time slot

  const Locations_ & locations_;       /// locations of observations
  const Variables geovars_;            /// Variables needed from model
  GetValues_ getvals_;                 /// GetValues used to fill in GeoVaLs
  std::unique_ptr<GeoVaLs_> geovals_;  /// GeoVaLs that are filled in
  std::vector<size_t> sizes_;          /// Sizes (e.g. number of vertical levels)
                                       /// for all Variables in GeoVaLs
  bool initialized_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GetValuePost<MODEL, OBS>::GetValuePost(const eckit::Configuration & conf, const Geometry_ & geom,
                                       const util::DateTime & bgn, const util::DateTime & end,
                                       const Locations_ & locations, const Variables & vars)
  : winbgn_(bgn), winend_(end), hslot_(), locations_(locations), geovars_(vars),
    getvals_(geom, locations_, conf), geovals_(), sizes_(geom.variableSizes(geovars_)),
    initialized_(false)
{
  Log::trace() << "GetValuePost::GetValuePost" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuePost<MODEL, OBS>::initialize(const util::Duration & tstep) {
  Log::trace() << "GetValuePost::doInitialize start" << std::endl;
  hslot_ = tstep/2;
  geovals_.reset(new GeoVaLs_(locations_, geovars_, sizes_));
  initialized_ = true;
  Log::trace() << "GetValuePost::doInitialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValuePost<MODEL, OBS>::process(const State_ & xx) {
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
