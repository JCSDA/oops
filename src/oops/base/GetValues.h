/*
 * (C) Copyright 2020-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredInterpolator.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/LocalInterpolator.h"
#include "oops/interface/Locations.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Timer.h"

namespace oops {


/// \brief Detect if \c MODEL defines \c LocalInterpolator.
///
/// If it does, HasInterpolator_ will be defined as std::true_type, otherwise as std::false_type.
///
/// \note Here's how this works. Note first that the primary template is followed by its
/// specialization. Suppose `HasInterpolator_ `is instantiated with a single parameter `MODEL`.
/// If `MODEL` doesn't define a type called `LocalInterpolator`, substitution of `MODEL` into the
/// specialisation will fail, so the primary template will be used. If `MODEL` defines
/// `LocalInterpolator`, the substitution of `MODEL` into `std::enable_if<...>` will succeed and
/// `std::enable_if<...>::type` will be defined as `void`.
/// This is identical to the value of the second parameter in the primary template. With both
/// parameters of the specialization matching that of the primary template, the compiler will
/// choose the specialization over the primary template.

// Primary template
template<typename MODEL, typename = void>
struct HasInterpolator_ : std::false_type {};

// Specialization
template<typename MODEL>
struct HasInterpolator_<MODEL,
       typename std::enable_if<std::is_convertible<typename MODEL::LocalInterpolator*,
                                                   void*>::value>::type>
    : std::true_type {};

// -----------------------------------------------------------------------------

template<typename MODEL, bool THasInterpolator_>
struct TModelInterpolator_IfAvailableElseGenericInterpolator;

template<typename MODEL>
struct TModelInterpolator_IfAvailableElseGenericInterpolator<MODEL, false> {
  typedef UnstructuredInterpolator<MODEL> type;
};

template<typename MODEL>
struct TModelInterpolator_IfAvailableElseGenericInterpolator<MODEL, true> {
  typedef LocalInterpolator<MODEL> type;
};

/// \brief Resolved to \c oops::LocalInterpolator<MODEL> (wrapper for \c MODEL::LocalInterpolator)
/// if \c MODEL defines \c LocalInterpolator; otherwise resolved to
/// \c oops::UnstructuredInterpolator<MODEL>
template<typename MODEL>
using TModelInterpolator_IfAvailableElseGenericInterpolator_t =
typename TModelInterpolator_IfAvailableElseGenericInterpolator<MODEL,
                                      HasInterpolator_<MODEL>::value>::type;

// -----------------------------------------------------------------------------

/// \brief Fills GeoVaLs with requested variables at obs locations during model run

template <typename MODEL, typename OBS>
class GetValues {
  typedef Geometry<MODEL>           Geometry_;
  typedef GeoVaLs<OBS>              GeoVaLs_;
  typedef Increment<MODEL>          Increment_;
  typedef TModelInterpolator_IfAvailableElseGenericInterpolator_t<MODEL> LocalInterp_;
  typedef Locations<OBS>            Locations_;
  typedef State<MODEL>              State_;

 public:
  GetValues(const eckit::Configuration &, const Geometry_ &,
            const util::DateTime &, const util::DateTime &,
            const Locations_ &, const Variables &, const Variables & varl = Variables());

/// Nonlinear
  void initialize(const util::Duration &);
  void process(const State_ &);
  void finalize();
  std::unique_ptr<GeoVaLs_> releaseGeoVaLs();

/// TL
  void initializeTL(const util::Duration &);
  void processTL(const Increment_ &);
  void finalizeTL();
  std::unique_ptr<GeoVaLs_> releaseGeoVaLsTL();

/// AD
  void releaseGeoVaLsAD(std::unique_ptr<GeoVaLs_> &);
  void initializeAD();
  void processAD(Increment_ &);
  void finalizeAD(const util::Duration &);

/// Variables that will be required from the State and Increment
  const Variables & linearVariables() const {return linvars_;}
  const Variables & requiredVariables() const {return geovars_;}

 private:
  util::DateTime winbgn_;   /// Begining of assimilation window
  util::DateTime winend_;   /// End of assimilation window
  util::Duration hslot_;    /// Half time slot

  const Locations_ & locations_;       /// locations of observations
  const Variables geovars_;            /// Variables needed from model
  std::vector<size_t> varsizes_;       /// Sizes (e.g. number of vertical levels)
                                       /// for all Variables in GeoVaLs
  std::unique_ptr<GeoVaLs_> geovals_;  /// GeoVaLs that are filled in

  const Variables linvars_;
  const std::vector<size_t> linsizes_;
  std::unique_ptr<GeoVaLs_> gvalstl_;        /// GeoVaLs for TL
  std::unique_ptr<const GeoVaLs_> gvalsad_;  /// Input GeoVaLs for adjoint forcing
  eckit::LocalConfiguration interpConf_;
  const eckit::mpi::Comm & comm_;
  const size_t ntasks_;
  std::vector<std::unique_ptr<LocalInterp_>> interp_;
  std::vector<std::vector<size_t>> myobs_index_by_task_;
  std::vector<std::vector<util::DateTime>> obs_times_by_task_;
  std::vector<std::vector<double>> locinterp_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GetValues<MODEL, OBS>::GetValues(const eckit::Configuration & conf, const Geometry_ & geom,
                                 const util::DateTime & bgn, const util::DateTime & end,
                                 const Locations_ & locs,
                                 const Variables & vars, const Variables & varl)
  : winbgn_(bgn), winend_(end), hslot_(), locations_(locs),
    geovars_(vars), varsizes_(geom.variableSizes(geovars_)), geovals_(),
    linvars_(varl), linsizes_(geom.variableSizes(linvars_)), gvalstl_(), gvalsad_(),
    interpConf_(conf), comm_(geom.getComm()), ntasks_(comm_.size()), interp_(ntasks_),
    myobs_index_by_task_(ntasks_), obs_times_by_task_(ntasks_), locinterp_()
{
  Log::trace() << "GetValues::GetValues start" << std::endl;
  util::Timer timer("oops::GetValues", "GetValues");

// Local obs coordinates
  std::vector<double> obslats = locations_.latitudes();
  std::vector<double> obslons = locations_.longitudes();
  std::vector<util::DateTime> obstimes = locations_.times();

// Exchange obs locations
  std::vector<std::vector<double>> myobs_locs_by_task(ntasks_);
  for (size_t jobs = 0; jobs < obstimes.size(); ++jobs) {
    const size_t itask = geom.closestTask(obslats[jobs], obslons[jobs]);
    myobs_index_by_task_[itask].push_back(jobs);
    myobs_locs_by_task[itask].push_back(obslats[jobs]);
    myobs_locs_by_task[itask].push_back(obslons[jobs]);
    obstimes[jobs].serialize(myobs_locs_by_task[itask]);
  }

  std::vector<std::vector<double>> mylocs_by_task(ntasks_);
  comm_.allToAll(myobs_locs_by_task, mylocs_by_task);

// Setup interpolators
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    // The 4 below is because each loc holds lat + lon + 2 datetime ints
    const size_t nobs = mylocs_by_task[jtask].size() / 4;
    std::vector<double> lats(nobs);
    std::vector<double> lons(nobs);
    obs_times_by_task_[jtask].resize(nobs);
    size_t ii = 0;
    for (size_t jobs = 0; jobs < nobs; ++jobs) {
      lats[jobs] = mylocs_by_task[jtask][ii];
      lons[jobs] = mylocs_by_task[jtask][ii + 1];
      ii += 2;
      obs_times_by_task_[jtask][jobs].deserialize(mylocs_by_task[jtask], ii);
    }
    ASSERT(mylocs_by_task[jtask].size() == ii);
    interp_[jtask] = std::make_unique<LocalInterp_>(interpConf_, geom, lats, lons);
  }

  Log::trace() << "GetValues::GetValues done" << std::endl;
}

// -----------------------------------------------------------------------------
//  Forward methods (called from nonlinear run)
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::initialize(const util::Duration & tstep) {
  Log::trace() << "GetValues::initialize start" << std::endl;
  double missing = 0.0;
  missing = util::missingValue(missing);
  size_t nvars = 0;
  for (size_t jj = 0; jj < varsizes_.size(); ++jj) nvars += varsizes_[jj];
  ASSERT(locinterp_.empty());
  locinterp_.resize(ntasks_);
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    locinterp_[jtask].resize(obs_times_by_task_[jtask].size() * nvars, missing);
  }
  hslot_ = tstep / 2;
  Log::trace() << "GetValues::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::process(const State_ & xx) {
  Log::trace() << "GetValues::process start" << std::endl;
  util::Timer timer("oops::GetValues", "process");

  util::DateTime t1 = std::max(xx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(xx.validTime()+hslot_, winend_);

  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
//  Mask obs outside time slot
    std::vector<bool> mask(obs_times_by_task_[jtask].size());
    for (size_t jobs = 0; jobs < obs_times_by_task_[jtask].size(); ++jobs) {
      mask[jobs] = obs_times_by_task_[jtask][jobs] > t1 && obs_times_by_task_[jtask][jobs] <= t2;
    }

//  Local interpolation
    interp_[jtask]->apply(geovars_, xx, mask, locinterp_[jtask]);
  }

  Log::trace() << "GetValues::process done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::finalize() {
  Log::trace() << "GetValues::finalize start" << std::endl;
  util::Timer timer("oops::GetValues", "finalize");

// Send/receive interpolated values
  std::vector<std::vector<double>> recvinterp(ntasks_);
  comm_.allToAll(locinterp_, recvinterp);
  locinterp_.clear();

// Store output in GeoVaLs
  size_t nvars = 0;
  for (size_t jj = 0; jj < varsizes_.size(); ++jj) nvars += varsizes_[jj];
  ASSERT(!geovals_);
  geovals_.reset(new GeoVaLs_(locations_, geovars_, varsizes_));
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    const size_t nvals = myobs_index_by_task_[jtask].size() * nvars;
    ASSERT(recvinterp[jtask].size() == nvals);
    geovals_->fill(myobs_index_by_task_[jtask], recvinterp[jtask]);
  }

  Log::trace() << "GetValues::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<GeoVaLs<OBS>> GetValues<MODEL, OBS>::releaseGeoVaLs() {
  Log::trace() << "GetValues::releaseGeoVaLs" << std::endl;
  ASSERT(geovals_);
// Release ownership of GeoVaLs
  return std::move(geovals_);
}

// -----------------------------------------------------------------------------
//  TL methods
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::initializeTL(const util::Duration & tstep) {
  Log::trace() << "GetValues::initializeTL start" << std::endl;
  double missing = 0.0;
  missing = util::missingValue(missing);
  size_t nvars = 0;
  for (size_t jj = 0; jj < linsizes_.size(); ++jj) nvars += linsizes_[jj];
  ASSERT(locinterp_.empty());
  locinterp_.resize(ntasks_);
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    locinterp_[jtask].resize(obs_times_by_task_[jtask].size() * nvars, missing);
  }
  hslot_ = tstep/2;
  Log::trace() << "GetValues::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::processTL(const Increment_ & dx) {
  Log::trace() << "GetValues::processTL start" << std::endl;
  util::Timer timer("oops::GetValues", "processTL");

  util::DateTime t1 = std::max(dx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(dx.validTime()+hslot_, winend_);

  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
//  Mask obs outside time slot
    std::vector<bool> mask(obs_times_by_task_[jtask].size());
    for (size_t jobs = 0; jobs < obs_times_by_task_[jtask].size(); ++jobs) {
      mask[jobs] = obs_times_by_task_[jtask][jobs] > t1 && obs_times_by_task_[jtask][jobs] <= t2;
    }

//  Local interpolation
    interp_[jtask]->apply(linvars_, dx, mask, locinterp_[jtask]);
  }

  Log::trace() << "GetValues::processTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::finalizeTL() {
  Log::trace() << "GetValues::finalizeTL start" << std::endl;
  util::Timer timer("oops::GetValues", "finalizeTL");

// Send/receive interpolated values
  std::vector<std::vector<double>> recvinterp(ntasks_);
  comm_.allToAll(locinterp_, recvinterp);
  locinterp_.clear();

// Store output in GeoVaLs
  size_t nvars = 0;
  for (size_t jj = 0; jj < linsizes_.size(); ++jj) nvars += linsizes_[jj];
  ASSERT(!gvalstl_);
  gvalstl_.reset(new GeoVaLs_(locations_, linvars_, linsizes_));
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    const size_t nvals = myobs_index_by_task_[jtask].size() * nvars;
    ASSERT(recvinterp[jtask].size() == nvals);
    gvalstl_->fill(myobs_index_by_task_[jtask], recvinterp[jtask]);
  }

  Log::trace() << "GetValues::processTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<GeoVaLs<OBS>> GetValues<MODEL, OBS>::releaseGeoVaLsTL() {
  Log::trace() << "GetValues::finalizeTL" << std::endl;
  ASSERT(gvalstl_);
  // Release ownership of GeoVaLs
  return std::move(gvalstl_);
}

// -----------------------------------------------------------------------------
//  AD methods
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::initializeAD() {
  Log::trace() << "GetValues::initializeAD start" << std::endl;
  locinterp_.clear();
  Log::trace() << "GetValues::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::processAD(Increment_ & dx) {
  Log::trace() << "GetValues::processAD start" << std::endl;
  util::Timer timer("oops::GetValues", "processAD");

  util::DateTime t1 = std::max(dx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(dx.validTime()+hslot_, winend_);

  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
//  Mask obs outside time slot
    std::vector<bool> mask(obs_times_by_task_[jtask].size());
    for (size_t jobs = 0; jobs < obs_times_by_task_[jtask].size(); ++jobs) {
      mask[jobs] = obs_times_by_task_[jtask][jobs] > t1 && obs_times_by_task_[jtask][jobs] <= t2;
    }

//  (Adjoint of) Local interpolation
    interp_[jtask]->applyAD(linvars_, dx, mask, locinterp_[jtask]);
  }

  Log::trace() << "GetValues::processAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::finalizeAD(const util::Duration & tstep) {
  Log::trace() << "GetValues::finalizeAD start" << std::endl;
  util::Timer timer("oops::GetValues", "finalizeAD");

  double missing = 0.0;
  missing = util::missingValue(missing);
  size_t nvars = 0;
  for (size_t jj = 0; jj < linsizes_.size(); ++jj) nvars += linsizes_[jj];
  ASSERT(locinterp_.empty());
  locinterp_.resize(ntasks_);
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    locinterp_[jtask].resize(obs_times_by_task_[jtask].size() * nvars, missing);
  }
  hslot_ = tstep/2;

  ASSERT(gvalsad_);

// (Adjoint of) Store output in GeoVaLs
  std::vector<std::vector<double>> recvinterp(ntasks_);
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    const size_t nvals = myobs_index_by_task_[jtask].size() * nvars;
    recvinterp[jtask].resize(nvals);
    gvalsad_->fillAD(myobs_index_by_task_[jtask], recvinterp[jtask]);
  }

// (Adjoint of) Send/receive interpolated values
  comm_.allToAll(recvinterp, locinterp_);

  gvalsad_.reset();

  Log::trace() << "GetValues::processAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::releaseGeoVaLsAD(std::unique_ptr<GeoVaLs_> & geovals) {
// Take ownership of GeoVaLs
  ASSERT(!gvalsad_);
  gvalsad_ = std::move(geovals);
  Log::trace() << "GetValues::releaseGeoVaLsAD" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

