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

#include "atlas/field.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/LocalUnstructuredInterpolator.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/LocalInterpolator.h"
#include "oops/interface/Locations.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
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
  typedef LocalUnstructuredInterpolator<MODEL> type;
};

template<typename MODEL>
struct TModelInterpolator_IfAvailableElseGenericInterpolator<MODEL, true> {
  typedef LocalInterpolator<MODEL> type;
};

/// \brief Resolved to \c oops::LocalInterpolator<MODEL> (wrapper for \c MODEL::LocalInterpolator)
/// if \c MODEL defines \c LocalInterpolator; otherwise resolved to
/// \c oops::LocalUnstructuredInterpolator<MODEL>
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

// *********************************************************************************
// *** TEMPORARY INTERFACES: Nonlinear, trajectory and TL are the same code with ***
// *** TEMPORARY INTERFACES: different variables, will be collapsed together     ***
// *********************************************************************************

/// Nonlinear
  void initialize(const util::Duration &);
  void process(const State_ &);
  std::unique_ptr<GeoVaLs_> releaseGeoVaLs();

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

/// Variables that will be required from the State and Increment
  const Variables & linearVariables() const {return linvars_;}
  const Variables & requiredVariables() const {return geovars_;}

 private:
  void setTree(const Geometry_ &);
  int closestTask(const double &, const double &) const;

  util::DateTime winbgn_;        /// Begining of assimilation window
  util::DateTime winend_;        /// End of assimilation window
  util::Duration hslot_;         /// Half time slot
  bool duplicatedDistribution_;  /// Set true if your model duplicates fields
                                 /// (locations will not be redistributed)

  const Locations_ & locations_;       /// locations of observations
  const Variables geovars_;            /// Variables needed from model
  std::vector<size_t> varsizes_;       /// Sizes (e.g. number of vertical levels)
                                       /// for all Variables in GeoVaLs
  std::unique_ptr<GeoVaLs_> geovals_;  /// GeoVaLs that are filled in

  const Variables linvars_;
  const std::vector<size_t> linsizes_;
  std::unique_ptr<GeoVaLs_> gvalstl_;        /// GeoVaLs for TL
  std::unique_ptr<const GeoVaLs_> gvalsad_;  /// Input GeoVaLs for adjoint forcing

  const eckit::mpi::Comm & comm_;
  const atlas::Geometry earth_;
  atlas::util::IndexKDTree globalTree_;
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
    comm_(geom.getComm()), earth_(atlas::util::Earth::radius()), globalTree_(earth_),
    duplicatedDistribution_(conf.getBool("duplicated distribution", false))
{
  Log::trace() << "GetValues::GetValues start" << std::endl;
  util::Timer timer("oops::GetValues", "GetValues");
  if (!duplicatedDistribution_) this->setTree(geom);
  Log::trace() << "GetValues::GetValues done" << std::endl;
}

// -----------------------------------------------------------------------------
//  Forward methods (called from nonlinear run)
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::initialize(const util::Duration & tstep) {
  Log::trace() << "GetValues::initialize start" << std::endl;
  hslot_ = tstep/2;
  geovals_.reset(new GeoVaLs_(locations_, geovars_, varsizes_));
  Log::trace() << "GetValues::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::process(const State_ & xx) {
  Log::trace() << "GetValues::process start" << std::endl;
  util::Timer timer("oops::GetValues", "process");

  eckit::LocalConfiguration conf;
  ASSERT(geovals_);
  size_t nvars = 0;
  for (size_t jj = 0; jj < varsizes_.size(); ++jj) nvars += varsizes_[jj];
  util::DateTime t1 = std::max(xx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(xx.validTime()+hslot_, winend_);
  const size_t ntasks = comm_.size();

// Get local obs coordinates
  std::vector<double> obslats;
  std::vector<double> obslons;
  std::vector<size_t> obsindx;
  locations_.localCoords(t1, t2, obslats, obslons, obsindx);

// Exchange obs locations
  std::vector<std::vector<size_t>> myobs_index_by_task(ntasks);
  std::vector<std::vector<double>> myobs_latlon_by_task(ntasks);
  for (size_t jobs = 0; jobs < obsindx.size(); ++jobs) {
    const size_t itask = duplicatedDistribution_ ? comm_.rank() :
                                                   this->closestTask(obslats[jobs], obslons[jobs]);
    myobs_index_by_task[itask].push_back(obsindx[jobs]);
    myobs_latlon_by_task[itask].push_back(obslats[jobs]);
    myobs_latlon_by_task[itask].push_back(obslons[jobs]);
  }
  std::vector<std::vector<double>> mylocs_latlon_by_task(ntasks);
  comm_.allToAll(myobs_latlon_by_task, mylocs_latlon_by_task);

// Interpolate
  std::vector<std::vector<double>> locinterp(ntasks);
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    LocalInterp_ interp(conf, xx.geometry(), mylocs_latlon_by_task[jtask]);
    interp.apply(geovars_, xx, locinterp[jtask]);
  }

// Send/receive interpolated values
  std::vector<std::vector<double>> recvinterp(ntasks);
  comm_.allToAll(locinterp, recvinterp);

// Store output in GeoVaLs
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    const size_t nvals = myobs_index_by_task[jtask].size() * nvars;
    ASSERT(recvinterp[jtask].size() == nvals);
    geovals_->fill(myobs_index_by_task[jtask], recvinterp[jtask]);
  }

  Log::trace() << "GetValues::process done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<GeoVaLs<OBS>> GetValues<MODEL, OBS>::releaseGeoVaLs() {
  Log::trace() << "GetValues::releaseGeoVaLs" << std::endl;
// Release ownership of GeoVaLs
  return std::move(geovals_);
}

// -----------------------------------------------------------------------------
//  Trajectory methods
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::initializeTraj(const util::Duration & tstep) {
  Log::trace() << "GetValues::initializeTraj start" << std::endl;
  hslot_ = tstep/2;
  geovals_.reset(new GeoVaLs_(locations_, geovars_, varsizes_));
  Log::trace() << "GetValues::initializeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::processTraj(const State_ & xx) {
  Log::trace() << "GetValues::processTraj start" << std::endl;
  util::Timer timer("oops::GetValues", "processTraj");

  eckit::LocalConfiguration conf;
  ASSERT(geovals_);
  size_t nvars = 0;
  for (size_t jj = 0; jj < varsizes_.size(); ++jj) nvars += varsizes_[jj];
  util::DateTime t1 = std::max(xx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(xx.validTime()+hslot_, winend_);
  const eckit::mpi::Comm & comm = xx.geometry().getComm();
  const size_t ntasks = comm.size();

// Get local obs coordinates
  std::vector<double> obslats;
  std::vector<double> obslons;
  std::vector<size_t> obsindx;
  locations_.localCoords(t1, t2, obslats, obslons, obsindx);

// Exchange obs locations
  std::vector<std::vector<size_t>> myobs_index_by_task(ntasks);
  std::vector<std::vector<double>> myobs_latlon_by_task(ntasks);
  for (size_t jobs = 0; jobs < obsindx.size(); ++jobs) {
    const size_t itask = this->closestTask(obslats[jobs], obslons[jobs]);
    myobs_index_by_task[itask].push_back(obsindx[jobs]);
    myobs_latlon_by_task[itask].push_back(obslats[jobs]);
    myobs_latlon_by_task[itask].push_back(obslons[jobs]);
  }
  std::vector<std::vector<double>> mylocs_latlon_by_task(ntasks);
  comm.allToAll(myobs_latlon_by_task, mylocs_latlon_by_task);

// Interpolate
  std::vector<std::vector<double>> locinterp(ntasks);
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    LocalInterp_ interp(conf, xx.geometry(), mylocs_latlon_by_task[jtask]);
    interp.apply(geovars_, xx, locinterp[jtask]);
  }

// Send/receive interpolated values
  std::vector<std::vector<double>> recvinterp(ntasks);
  comm.allToAll(locinterp, recvinterp);

// Store output in GeoVaLs
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    const size_t nvals = myobs_index_by_task[jtask].size() * nvars;
    ASSERT(recvinterp[jtask].size() == nvals);
    geovals_->fill(myobs_index_by_task[jtask], recvinterp[jtask]);
  }

  Log::trace() << "GetValues::processTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<GeoVaLs<OBS>> GetValues<MODEL, OBS>::finalizeTraj() {
  Log::trace() << "GetValues::finalizeTraj" << std::endl;
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
  hslot_ = tstep/2;
  gvalstl_.reset(new GeoVaLs_(locations_, linvars_, linsizes_));
  Log::trace() << "GetValues::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::processTL(const Increment_ & dx) {
  Log::trace() << "GetValues::processTL start" << std::endl;
  util::Timer timer("oops::GetValues", "processTL");

  eckit::LocalConfiguration conf;
  ASSERT(gvalstl_);
  size_t nvars = 0;
  for (size_t jj = 0; jj < linsizes_.size(); ++jj) nvars += linsizes_[jj];
  util::DateTime t1 = std::max(dx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(dx.validTime()+hslot_, winend_);
  const eckit::mpi::Comm & comm = dx.geometry().getComm();
  const size_t ntasks = comm.size();

// Get local obs coordinates
  std::vector<double> obslats;
  std::vector<double> obslons;
  std::vector<size_t> obsindx;
  locations_.localCoords(t1, t2, obslats, obslons, obsindx);

// Exchange obs locations
  std::vector<std::vector<size_t>> myobs_index_by_task(ntasks);
  std::vector<std::vector<double>> myobs_latlon_by_task(ntasks);
  for (size_t jobs = 0; jobs < obsindx.size(); ++jobs) {
    const size_t itask = this->closestTask(obslats[jobs], obslons[jobs]);
    myobs_index_by_task[itask].push_back(obsindx[jobs]);
    myobs_latlon_by_task[itask].push_back(obslats[jobs]);
    myobs_latlon_by_task[itask].push_back(obslons[jobs]);
  }
  std::vector<std::vector<double>> mylocs_latlon_by_task(ntasks);
  comm.allToAll(myobs_latlon_by_task, mylocs_latlon_by_task);

// Interpolate
  std::vector<std::vector<double>> locinterp(ntasks);
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    LocalInterp_ interp(conf, dx.geometry(), mylocs_latlon_by_task[jtask]);
    interp.apply(linvars_, dx, locinterp[jtask]);
  }

// Send/receive interpolated values
  std::vector<std::vector<double>> recvinterp(ntasks);
  comm.allToAll(locinterp, recvinterp);

// Store output in GeoVaLs
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    const size_t nvals = myobs_index_by_task[jtask].size() * nvars;
    ASSERT(recvinterp[jtask].size() == nvals);
    gvalstl_->fill(myobs_index_by_task[jtask], recvinterp[jtask]);
  }

  Log::trace() << "GetValues::processTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<GeoVaLs<OBS>> GetValues<MODEL, OBS>::finalizeTL() {
  Log::trace() << "GetValues::finalizeTL" << std::endl;
  ASSERT(gvalstl_);
  // Release ownership of GeoVaLs
  return std::move(gvalstl_);
}

// -----------------------------------------------------------------------------
//  AD methods
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::setAD(std::unique_ptr<GeoVaLs_> & geovals) {
// Take ownership of GeoVaLs
  gvalsad_ = std::move(geovals);
  Log::trace() << "GetValues::setAD" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::initializeAD(const util::Duration & tstep) {
  hslot_ = tstep/2;
  Log::trace() << "GetValues::initializeAD" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::processAD(Increment_ & dx) {
  Log::trace() << "GetValues::processAD start" << std::endl;
  util::Timer timer("oops::GetValues", "processAD");

  eckit::LocalConfiguration conf;
  ASSERT(gvalsad_);
  size_t nvars = 0;
  for (size_t jj = 0; jj < linsizes_.size(); ++jj) nvars += linsizes_[jj];
  util::DateTime t1 = std::max(dx.validTime()-hslot_, winbgn_);
  util::DateTime t2 = std::min(dx.validTime()+hslot_, winend_);
  const eckit::mpi::Comm & comm = dx.geometry().getComm();
  const size_t ntasks = comm.size();

// Get local obs coordinates
  std::vector<double> obslats;
  std::vector<double> obslons;
  std::vector<size_t> obsindx;
  locations_.localCoords(t1, t2, obslats, obslons, obsindx);

// Exchange obs locations
  std::vector<std::vector<size_t>> myobs_index_by_task(ntasks);
  std::vector<std::vector<double>> myobs_latlon_by_task(ntasks);
  for (size_t jobs = 0; jobs < obsindx.size(); ++jobs) {
    const size_t itask = this->closestTask(obslats[jobs], obslons[jobs]);
    myobs_index_by_task[itask].push_back(obsindx[jobs]);
    myobs_latlon_by_task[itask].push_back(obslats[jobs]);
    myobs_latlon_by_task[itask].push_back(obslons[jobs]);
  }
  std::vector<std::vector<double>> mylocs_latlon_by_task(ntasks);
  comm.allToAll(myobs_latlon_by_task, mylocs_latlon_by_task);

// (Adjoint of) Store output in GeoVaLs
  std::vector<std::vector<double>> recvinterp(ntasks);
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    const size_t nvals = myobs_index_by_task[jtask].size() * nvars;
    recvinterp[jtask].resize(nvals);
    gvalsad_->fillAD(myobs_index_by_task[jtask], recvinterp[jtask]);
  }

// (Adjoint of) Send/receive interpolated values
  std::vector<std::vector<double>> locinterp(ntasks);
  comm.allToAll(recvinterp, locinterp);

// (Adjoint of) Interpolate
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    LocalInterp_ interp(conf, dx.geometry(), mylocs_latlon_by_task[jtask]);
    interp.applyAD(linvars_, dx, locinterp[jtask]);
  }

  Log::trace() << "GetValues::processAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::finalizeAD() {
  ASSERT(gvalsad_);
  gvalsad_.reset();
}

// -----------------------------------------------------------------------------
//  Private methods
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::setTree(const Geometry_ & grid) {
  const size_t ntasks = comm_.size();

// Local latitudes and longitudes
  std::vector<double> lats;
  std::vector<double> lons;
  grid.latlon(lats, lons, false);
  for (size_t jj = 0; jj < lons.size(); ++jj) {
    if (lons[jj] < 0.0) lons[jj] += 360.0;
  }

  const size_t sizel = lats.size();
  std::vector<double> latlon(2*sizel);
  for (size_t jj = 0; jj < sizel; ++jj) {
    latlon[2*jj]   = lats[jj];
    latlon[2*jj+1] = lons[jj];
  }

// Collect global grid lats and lons
  std::vector<size_t> sizes(ntasks);
  comm_.allGather(sizel, sizes.begin(), sizes.end());

  size_t sizeg = 0;
  for (size_t jtask = 0; jtask < ntasks; ++jtask) sizeg += sizes[jtask];

  std::vector<double> latlon_global(2*sizeg);
  mpi::allGatherv(comm_, latlon, latlon_global);

// Arrange coordinates and task index for kd-tree
  std::vector<double> latglo(sizeg), longlo(sizeg);
  std::vector<int> taskindx(sizeg);
  size_t jglo = 0;
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    for (size_t jj = 0; jj < sizes[jtask]; ++jj) {
      latglo[jglo] = latlon_global[2*jglo];
      longlo[jglo] = latlon_global[2*jglo+1];
      taskindx[jglo] = jtask;
      ++jglo;
    }
  }
  ASSERT(jglo == sizeg);

// Create global kd-tree
  globalTree_.build(longlo, latglo, taskindx);

  Log::info() << "GetValues: Global tree size = " << globalTree_.size()
              << ", footprint = " << globalTree_.footprint() << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
int GetValues<MODEL, OBS>::closestTask(const double & lat, const double & lon) const {
  atlas::PointLonLat obsloc(lon, lat);
  obsloc.normalise();
  const int itask = globalTree_.closestPoint(obsloc).payload();
  ASSERT(itask >= 0 && itask < comm_.size());
  return itask;
}

// -----------------------------------------------------------------------------

}  // namespace oops

