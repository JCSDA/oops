/*
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_GETVALUETLAD_H_
#define OOPS_BASE_GETVALUETLAD_H_

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/LinearGetValues.h"
#include "oops/interface/Locations.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// GetValueTLAD is its own class for now. Once the change of variables is moved out it
// should be a subclass of GetValuePost that adds setGeoVaLs, processTL and processAD
// The other methods are the same.

/// \brief TLAD of filling GeoVaLs with requested variables at requested locations during model run
template <typename MODEL, typename OBS>
class GetValueTLAD {
  typedef Geometry<MODEL>              Geometry_;
  typedef GeoVaLs<OBS>                 GeoVaLs_;
  typedef LinearGetValues<MODEL, OBS>  GetValues_;
  typedef Increment<MODEL>             Increment_;
  typedef Locations<OBS>               Locations_;
  typedef State<MODEL>                 State_;

 public:
/// \brief Saves Locations and Variables to be processed
  GetValueTLAD(const eckit::Configuration &, const Geometry_ &,
               const util::DateTime &, const util::DateTime &,
               const Locations_ &, const Variables &, const Variables &);

/// Linearization trajectory
  void initializeTraj(const util::Duration &);
  void processTraj(const State_ &);
  std::unique_ptr<GeoVaLs_> finalizeTraj();

/// TL
  void initializeTL(const util::Duration &);
  void processTL(const Increment_ &);
  std::unique_ptr<GeoVaLs_> finalizeTL();

/// AD
  void setAD(std::unique_ptr<GeoVaLs_> &);
  void initializeAD(const util::Duration &);
  void processAD(Increment_ &);
  void finalizeAD();

/// Variables that will be required from the State
  const Variables & requiredVariables() const {return geovars_;}
  const Variables & linearVariables() const {return linvars_;}

 private:
  util::DateTime winbgn_;   /// Begining of assimilation window
  util::DateTime winend_;   /// End of assimilation window
  util::Duration hslot_;    /// Half time slot

  const Locations_ & locations_;             /// locations of observations
  const Variables geovars_;                  /// Variables needed from model
  const std::vector<size_t> geovars_sizes_;  /// Sizes (e.g. number of vertical levels)
                                             /// for all geovars_ variables
  const Variables linvars_;                  /// Variables needed from linear model
  const std::vector<size_t> linvars_sizes_;  /// Sizes (e.g. number of vertical levels)
                                             /// for all linvars_ variables
  GetValues_ getvals_;                       /// GetValues used to fill in GeoVaLs
  std::unique_ptr<GeoVaLs_> geovals_;        /// GeoVaLs that are filled in
  std::unique_ptr<GeoVaLs_> gvalstl_;        /// GeoVaLs for TL
  std::unique_ptr<const GeoVaLs_> gvalsad_;  /// Input GeoVaLs for adjoint forcing
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GetValueTLAD<MODEL, OBS>::GetValueTLAD(const eckit::Configuration & conf, const Geometry_ & geom,
                                       const util::DateTime & bgn, const util::DateTime & end,
                                       const Locations_ & locations,
                                       const Variables & vars, const Variables & varl)
  : winbgn_(bgn), winend_(end), hslot_(), locations_(locations),
    geovars_(vars), geovars_sizes_(geom.variableSizes(geovars_)),
    linvars_(varl), linvars_sizes_(geom.variableSizes(linvars_)),
    getvals_(geom, locations_, conf), geovals_(), gvalstl_(), gvalsad_()
{
  Log::trace() << "GetValueTLAD::GetValueTLAD" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValueTLAD<MODEL, OBS>::initializeTraj(const util::Duration & tstep) {
  Log::trace() << "GetValueTLAD::initializeTraj start" << std::endl;
  hslot_ = tstep/2;
  geovals_.reset(new GeoVaLs_(locations_, geovars_, geovars_sizes_));
  Log::trace() << "GetValueTLAD::initializeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValueTLAD<MODEL, OBS>::processTraj(const State_ & xx) {
  Log::trace() << "GetValueTLAD::processTraj start" << std::endl;
  ASSERT(geovals_);
  util::DateTime t1 = std::max(xx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(xx.validTime()+hslot_, winend_);

// Get state variables at obs locations
  getvals_.setTrajectory(xx, t1, t2, *geovals_);
  Log::trace() << "GetValueTLAD::processTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<GeoVaLs<OBS>> GetValueTLAD<MODEL, OBS>::finalizeTraj() {
  Log::trace() << "GetValueTLAD::finalizeTraj" << std::endl;
  ASSERT(geovals_);
  // Release ownership of GeoVaLs
  return std::move(geovals_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValueTLAD<MODEL, OBS>::initializeTL(const util::Duration & tstep) {
  Log::trace() << "GetValueTLAD::initializeTL start" << std::endl;
  hslot_ = tstep/2;
  gvalstl_.reset(new GeoVaLs_(locations_, linvars_, linvars_sizes_));
  Log::trace() << "GetValueTLAD::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValueTLAD<MODEL, OBS>::processTL(const Increment_ & dx) {
  Log::trace() << "GetValueTLAD::processTL start" << std::endl;
  ASSERT(gvalstl_);
  util::DateTime t1 = std::max(dx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(dx.validTime()+hslot_, winend_);

// Get state variables at obs locations
  getvals_.fillGeoVaLsTL(dx, t1, t2, *gvalstl_);
  Log::trace() << "GetValueTLAD::processTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<GeoVaLs<OBS>> GetValueTLAD<MODEL, OBS>::finalizeTL() {
  Log::trace() << "GetValueTLAD::finalizeTL" << std::endl;
  ASSERT(gvalstl_);
  // Release ownership of GeoVaLs
  return std::move(gvalstl_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValueTLAD<MODEL, OBS>::setAD(std::unique_ptr<GeoVaLs_> & geovals) {
  // Take ownership of GeoVaLs
  gvalsad_ = std::move(geovals);
  Log::trace() << "GetValueTLAD::setAD" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValueTLAD<MODEL, OBS>::initializeAD(const util::Duration & tstep) {
  hslot_ = tstep/2;
  Log::trace() << "GetValueTLAD::initializeAD" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValueTLAD<MODEL, OBS>::processAD(Increment_ & dx) {
  Log::trace() << "GetValueTLAD::processAD start" << std::endl;
  ASSERT(gvalsad_);
  util::DateTime t1 = std::max(dx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(dx.validTime()+hslot_, winend_);

// Get state variables at obs locations
  getvals_.fillGeoVaLsAD(dx, t1, t2, *gvalsad_);

  Log::trace() << "GetValueTLAD::processAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValueTLAD<MODEL, OBS>::finalizeAD() {
  ASSERT(gvalsad_);
  gvalsad_.reset();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_GETVALUETLAD_H_
