/*
 * (C) Copyright 2020-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <Eigen/Core>

#include <algorithm>
#include <memory>
#include <numeric>
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
#include "oops/interface/SampledLocations.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Timer.h"
#include "oops/util/TimeWindow.h"

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
  typedef UnstructuredInterpolator type;
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

// PreProcessModel provides a specialization that handles models whose model-specific interpolator
// acts on atlas::FieldSets (there we call the atlas halo exchanges), and a default case that
// handles models whose interpolator acts on model States and Increments (no-op).
template <typename MODEL, typename = cpp17::void_t<>>
struct PreProcessModel {
  static void preProcessModelData(const oops::State<MODEL> & state) { return; }
  static void preProcessModelData(const oops::Increment<MODEL> & increment) { return; }
  static void preProcessModelDataAD(const oops::Increment<MODEL> & increment) { return; }
};
template <typename MODEL>
struct PreProcessModel<MODEL, cpp17::void_t<
           decltype(std::declval<typename MODEL::LocalInterpolator>().
               apply(std::declval<oops::Variables>(),
                 std::declval<atlas::FieldSet>(),
                 std::declval<std::vector<bool>>(),
                 std::declval<std::vector<double>&>()))>> {
  static void preProcessModelData(const oops::State<MODEL> & state) {
    const atlas::FieldSet & fset = state.fieldSet().fieldSet();
    fset.haloExchange();
  }
  static void preProcessModelData(const oops::Increment<MODEL> & increment) {
    const atlas::FieldSet & fset = increment.fieldSet().fieldSet();
    fset.haloExchange();
  }
  static void preProcessModelDataAD(const oops::Increment<MODEL> & increment) {
    const atlas::FieldSet & fset = increment.fieldSet().fieldSet();
    fset.adjointHaloExchange();
    fset.set_dirty();
  }
};

// Switch between two pre-process implementations depending on whether the MODEL is using a
// specific interpolator or is falling back to the generic UnstructuredInterpolator.
template <typename MODEL, bool HasInterpolator>
struct PreProcessModelOrGeneric;
// This `false` specialization handles the oops UnstructuredInterpolator which works on generic
// data in atlas::FielSets; in this case we call the atlas halo exchange.
template <typename MODEL>
struct PreProcessModelOrGeneric<MODEL, false> {
  static void preProcessModelData(const oops::State<MODEL> & state) {
    const atlas::FieldSet & fset = state.fieldSet().fieldSet();
    fset.haloExchange();
  }
  static void preProcessModelData(const oops::Increment<MODEL> & increment) {
    const atlas::FieldSet & fset = increment.fieldSet().fieldSet();
    fset.haloExchange();
  }
  static void preProcessModelDataAD(const oops::Increment<MODEL> & increment) {
    const atlas::FieldSet & fset = increment.fieldSet().fieldSet();
    fset.adjointHaloExchange();
    fset.set_dirty();
  }
};
// This `true` specialization handles model-specific interpolators; we forward to another level of
// helper to switch between model-specific interpolators acting on atlas vs on model-specific data.
template <typename MODEL>
struct PreProcessModelOrGeneric<MODEL, true> {
  static void preProcessModelData(const oops::State<MODEL> & state) {
    PreProcessModel<MODEL>::preProcessModelData(state);
  }
  static void preProcessModelData(const oops::Increment<MODEL> & increment) {
    PreProcessModel<MODEL>::preProcessModelData(increment);
  }
  static void preProcessModelDataAD(const oops::Increment<MODEL> & increment) {
    PreProcessModel<MODEL>::preProcessModelDataAD(increment);
  }
};

// PreProcessHelper selects whether to preprocess the model State or Increment before passing it to
// the interpolator. This allows for calling an atlas haloExchange only in the case where the
// interpolator has an atlas-based apply method. See the similar logic within the LocalInterpolator
// itself; we duplicate this logic here because we want to keep the LocalInterpolator free of MPI.
template <typename MODEL>
struct PreProcessHelper {
  static void preProcessModelData(const oops::State<MODEL> & state) {
    PreProcessModelOrGeneric<MODEL, HasInterpolator_<MODEL>::value>::preProcessModelData(state);
  }
  static void preProcessModelData(const oops::Increment<MODEL> & increment) {
    PreProcessModelOrGeneric<MODEL, HasInterpolator_<MODEL>::value>::preProcessModelData(increment);
  }
  static void preProcessModelDataAD(const oops::Increment<MODEL> & increment) {
    PreProcessModelOrGeneric<MODEL, HasInterpolator_<MODEL>::value>::preProcessModelDataAD(
        increment);
  }
};

// -----------------------------------------------------------------------------

/// \brief Fills GeoVaLs with requested variables at obs locations during model run

template <typename MODEL, typename OBS>
class GetValues : private util::ObjectCounter<GetValues<MODEL, OBS> > {
  typedef Geometry<MODEL>           Geometry_;
  typedef GeoVaLs<OBS>              GeoVaLs_;
  typedef Increment<MODEL>          Increment_;
  typedef TModelInterpolator_IfAvailableElseGenericInterpolator_t<MODEL> LocalInterp_;
  typedef SampledLocations<OBS>     SampledLocations_;
  typedef State<MODEL>              State_;

 public:
  static const std::string classname() {return "oops::GetValues";}

  GetValues(const eckit::Configuration &, const Geometry_ &,
            const util::TimeWindow &,
            const SampledLocations_ &,
            const Variables &, const Variables & varl = Variables());

/// Nonlinear
  void initialize(const util::Duration &);
  void process(const State_ &);
  void finalize();
  void fillGeoVaLs(GeoVaLs_ &);

/// TL
  void initializeTL(const util::Duration &);
  void processTL(const Increment_ &);
  void finalizeTL();
  void fillGeoVaLsTL(GeoVaLs_ &);

/// AD
  void fillGeoVaLsAD(const GeoVaLs_ &);
  void initializeAD();
  void processAD(Increment_ &);
  void finalizeAD(const util::Duration &);

/// Variables that will be required from the State and Increment
  const Variables & linearVariables() const {return linvars_;}
  const Variables & requiredVariables() const {return geovars_;}
  const bool & useMethodsTL() const {return geovalsTL_;}

 private:
/// time-interpolation helper: adds contribution from this time to running total
  void incInterpValues(const util::DateTime &, const std::vector<bool> &,
                       const size_t &, const std::vector<double> &);


  util::Duration hslot_;    /// Half time slot
  const util::TimeWindow timeWindow_;

  const Variables geovars_;  /// Variables needed from model
  size_t varsizes_;          /// Sizes (e.g. number of vertical levels)
                             /// for all Variables in GeoVaLs
  const Variables linvars_;
  size_t linsizes_;
  eckit::LocalConfiguration interpConf_;
  const eckit::mpi::Comm & comm_;
  const size_t ntasks_;
  std::vector<std::unique_ptr<LocalInterp_>> interp_;
  std::vector<std::vector<size_t>> myobs_index_by_task_;
  std::vector<std::vector<util::DateTime>> obs_times_by_task_;
  std::vector<std::vector<double>> locinterp_;
  std::vector<std::vector<double>> recvinterp_;
  std::vector<eckit::mpi::Request> send_req_;
  std::vector<eckit::mpi::Request> recv_req_;
  int tag_;
  const bool levelsTopDown_;           /// When true: Levels are in top down order.
  std::vector<size_t> geovarsSizes_;   /// number of levels for geovars_
  std::vector<size_t> linvarsSizes_;   /// number of levels for linvars_
  bool useLinearTimeInterpolation_;    /// set true for linear and false for
                                       /// nearest-neighbour time-
                                       /// interpolation (default false)
  bool doLinearTimeInterpolation_;     /// set true when linear time interpolation
                                       /// needs to be done for this run
  std::vector<size_t> recv_tasks_;
  bool geovalsTL_ = false;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GetValues<MODEL, OBS>::GetValues(const eckit::Configuration & conf, const Geometry_ & geom,
                                 const util::TimeWindow & timeWindow,
                                 const SampledLocations_ & locs,
                                 const Variables & vars, const Variables & varl)
  : timeWindow_(timeWindow),
    geovars_(vars), varsizes_(0), linvars_(varl), linsizes_(0),
    interpConf_(conf), comm_(geom.getComm()), ntasks_(comm_.size()), interp_(ntasks_),
    myobs_index_by_task_(ntasks_), obs_times_by_task_(ntasks_),
    locinterp_(), recvinterp_(), send_req_(), recv_req_(), tag_(789),
    levelsTopDown_(geom.levelsAreTopDown()), geovarsSizes_(geom.variableSizes(geovars_)),
    linvarsSizes_(geom.variableSizes(linvars_))
{
  Log::trace() << "GetValues::GetValues start" << std::endl;
  util::Timer timer("oops::GetValues", "GetValues");

// set the type of time-interpolation
  std::string value;
  if (conf.get("time interpolation", value)) {
    if (value == "linear") {
      useLinearTimeInterpolation_ = true;
    } else if (value == "nearest") {
      useLinearTimeInterpolation_ = false;
    } else {
      ABORT("GetValues::GetValues: time interpolation has an unsupported value.");
    }
  } else {
    useLinearTimeInterpolation_ = false;
  }

  tag_ += this->created();

  varsizes_ = std::accumulate(geovarsSizes_.begin(), geovarsSizes_.end(), 0);
  linsizes_ = std::accumulate(linvarsSizes_.begin(), linvarsSizes_.end(), 0);

// Interpolation paths (currently vertical columns) sampling the obs locations

  const std::vector<double> &obslats = locs.latitudes();
  const std::vector<double> &obslons = locs.longitudes();
  const std::vector<util::DateTime> &obstimes = locs.times();

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
  const double missing = util::missingValue<double>();
  ASSERT(locinterp_.empty());

  locinterp_.resize(ntasks_);
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    locinterp_[jtask].resize(obs_times_by_task_[jtask].size() * varsizes_, missing);
  }
  doLinearTimeInterpolation_ = useLinearTimeInterpolation_;
  // no need to do time interpolation if there is only one subwindow
  if (tstep == timeWindow_.length()) doLinearTimeInterpolation_ = false;
  hslot_ = doLinearTimeInterpolation_ ? tstep : tstep/2;

  Log::trace() << "GetValues::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::incInterpValues(
                    const util::DateTime & tCurrent, const std::vector<bool> & mask,
                    const size_t & jtask,
                    const std::vector<double> & tmplocinterp)
{
  Log::trace() << "GetValues::incInterpValues start" << std::endl;

  const double missing = util::missingValue<double>();

// Get the state previous and next times and time-step
  const util::DateTime tPrevious = tCurrent - hslot_;
  const util::DateTime tNext = tCurrent + hslot_;
  const double dt = static_cast<double>(hslot_.toSeconds());

// Compute and add time weighted contribution from the input interpolated values.
  double timeWeight = 0;
  const size_t nObs = obs_times_by_task_[jtask].size();
  for (size_t jp = 0; jp < nObs; ++jp) {
    size_t valuesIndex = jp;
    if (mask[jp]) {
      const util::DateTime & obCurrentTime = obs_times_by_task_[jtask][jp];
      const bool isCurrentTime = obCurrentTime == tCurrent;
      const bool isFirst = obCurrentTime > tCurrent;
      if (!isCurrentTime && isFirst) {
        timeWeight =
          static_cast<double>((tNext - obCurrentTime).toSeconds())/dt;
      } else if (!isCurrentTime && !isFirst) {
        timeWeight =
          static_cast<double>((obCurrentTime - tPrevious).toSeconds())/dt;
      }
      for (size_t jf = 0; jf < geovars_.size(); ++jf) {
        for (size_t jlev = 0; jlev < geovarsSizes_[jf]; ++jlev) {
          if (tmplocinterp[valuesIndex] == missing) {
            locinterp_[jtask][valuesIndex] = missing;
          } else if (isCurrentTime) {
            locinterp_[jtask][valuesIndex] = tmplocinterp[valuesIndex];
          } else if (isFirst) {
            locinterp_[jtask][valuesIndex] = tmplocinterp[valuesIndex]*timeWeight;
          } else if (locinterp_[jtask][valuesIndex] != missing) {
            // Don't linearly interpolate missing data
            locinterp_[jtask][valuesIndex] += tmplocinterp[valuesIndex]*timeWeight;
          }
          valuesIndex += nObs;
        }
      }
    }
  }
  Log::trace() << "GetValues::incInterpValues done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::process(const State_ & xx) {
  Log::trace() << "GetValues::process start" << std::endl;
  util::Timer timer("oops::GetValues", "process");

  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
//  Mask obs outside time slot
    const util::TimeWindow timeSubWindow =
      timeWindow_.createSubWindow(xx.validTime(), hslot_);
    const std::vector<bool> mask =
      timeSubWindow.createTimeMask(obs_times_by_task_[jtask]);

//  Local interpolation
    if (doLinearTimeInterpolation_) {
      std::vector<double> tmplocinterp(locinterp_[jtask].size(), 0);
      interp_[jtask]->apply(geovars_, xx, mask, tmplocinterp);
      incInterpValues(xx.validTime(), mask, jtask, tmplocinterp);
    } else {
      interp_[jtask]->apply(geovars_, xx, mask, locinterp_[jtask]);
    }
  }

  Log::trace() << "GetValues::process done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::finalize() {
  Log::trace() << "GetValues::finalize start" << std::endl;
  util::Timer timer("oops::GetValues", "finalize");

// Allocate receive buffers and non blocking receive of interpolated values
  ASSERT(recvinterp_.empty());
  ASSERT(recv_req_.empty());
  ASSERT(recv_tasks_.empty());
  recvinterp_.resize(ntasks_);
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    const size_t nrecv = myobs_index_by_task_[jtask].size() * varsizes_;
    if (nrecv > 0) {
      recvinterp_[jtask].resize(nrecv);
      recv_req_.push_back(comm_.iReceive(&recvinterp_[jtask][0], nrecv, jtask, tag_));
      recv_tasks_.push_back(jtask);
    }
  }

// Send values interpolated locally (non-blocking)
  ASSERT(send_req_.empty());
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    if (locinterp_[jtask].size() > 0) {
      send_req_.push_back(comm_.iSend(&locinterp_[jtask][0], locinterp_[jtask].size(),
                                      jtask, tag_));
    }
  }

// Add MPI barrier to work around an intel MPI deadlock on some AMD platforms;
// barrier is added every N'th obs type
#ifdef INTELMPI_DEADLOCK_GETVALUES_LIMIT
  if ( tag_ % INTELMPI_DEADLOCK_GETVALUES_LIMIT == 0 ) comm_.barrier();
#endif

  Log::trace() << "GetValues::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::fillGeoVaLs(GeoVaLs_ & geovals) {
  Log::trace() << "GetValues::fillGeoVaLs start" << std::endl;
  util::Timer timer("oops::GetValues", "fillGeoVaLs");

// Wait for received interpolated values and store in GeoVaLs
  ASSERT(recvinterp_.size() == ntasks_);
  for (size_t jreq = 0; jreq < recv_req_.size(); ++jreq) {
    int ireq = -1;
    eckit::mpi::Status rst = comm_.waitAny(recv_req_, ireq);
    ASSERT(rst.error() == 0);
    size_t itask = recv_tasks_[ireq];
    ASSERT(itask >=0 && (size_t)itask < ntasks_);

    ASSERT(recvinterp_[itask].size() == myobs_index_by_task_[itask].size() * varsizes_);

    // Create non-owning views ("maps") into the interpolation results.
    const Eigen::Map<const Eigen::VectorX<size_t>> indices(myobs_index_by_task_[itask].data(),
                                                           myobs_index_by_task_[itask].size());
    // Each column contains the values of a single variable at a single level and all locations
    // with indices 'indices'. The columns are ordered first by level and then by variable.
    const Eigen::Map<const Eigen::MatrixXd> values(recvinterp_[itask].data(),
                                                   myobs_index_by_task_[itask].size(), varsizes_);

    size_t colOffset = 0;
    for (size_t jvar = 0; jvar < geovars_.size(); ++jvar) {
      const size_t numLevels = geovarsSizes_[jvar];
      geovals.fill(geovars_[jvar], indices, values.middleCols(colOffset, numLevels),
                   this->levelsTopDown_);
      colOffset += numLevels;
    }
    ASSERT(colOffset == varsizes_);
  }
  recv_req_.clear();
  recv_tasks_.clear();
  recvinterp_.clear();

// Clean-up send buffers (after making sure data has been sent)
  for (size_t jreq = 0; jreq < send_req_.size(); ++jreq) {
    int itask = -1;
    eckit::mpi::Status sst = comm_.waitAny(send_req_, itask);
    ASSERT(sst.error() == 0);
  }
  send_req_.clear();
  locinterp_.clear();

  Log::trace() << "GetValues::fillGeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------
//  TL methods
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::initializeTL(const util::Duration & tstep) {
  Log::trace() << "GetValues::initializeTL start" << std::endl;
  const double missing = util::missingValue<double>();
  ASSERT(locinterp_.empty());
  locinterp_.resize(ntasks_);
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    locinterp_[jtask].resize(obs_times_by_task_[jtask].size() * linsizes_, missing);
  }
  hslot_ = tstep/2;
  geovalsTL_ = true;
  Log::trace() << "GetValues::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::processTL(const Increment_ & dx) {
  Log::trace() << "GetValues::processTL start" << std::endl;
  util::Timer timer("oops::GetValues", "processTL");

  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
//  Mask obs outside time slot
    const util::TimeWindow timeSubWindow =
      timeWindow_.createSubWindow(dx.validTime(), hslot_);
    const std::vector<bool> mask =
      timeSubWindow.createTimeMask(obs_times_by_task_[jtask]);

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

// Allocate receive buffers and non blocking receive of interpolated values
  ASSERT(recvinterp_.empty());
  ASSERT(recv_req_.empty());
  ASSERT(recv_tasks_.empty());
  recvinterp_.resize(ntasks_);
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    const size_t nrecv = myobs_index_by_task_[jtask].size() * linsizes_;
    recvinterp_[jtask].resize(nrecv);
    if (nrecv > 0) {
      recv_req_.push_back(comm_.iReceive(&recvinterp_[jtask][0], nrecv, jtask, tag_));
      recv_tasks_.push_back(jtask);
    }
  }

// Send values interpolated locally (non-blocking)
  ASSERT(send_req_.empty());
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    if (locinterp_[jtask].size() > 0) {
      send_req_.push_back(comm_.iSend(&locinterp_[jtask][0], locinterp_[jtask].size(),
                                      jtask, tag_));
    }
  }

// Add MPI barrier to work around an intel MPI deadlock on some AMD platforms;
// barrier is added every N'th obs type
#ifdef INTELMPI_DEADLOCK_GETVALUES_LIMIT
  if ( tag_ % INTELMPI_DEADLOCK_GETVALUES_LIMIT == 0 ) comm_.barrier();
#endif

  Log::trace() << "GetValues::finalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::fillGeoVaLsTL(GeoVaLs_ & geovals) {
  Log::trace() << "GetValues::fillGeoVaLsTL start" << std::endl;
  util::Timer timer("oops::GetValues", "fillGeoVaLsTL");

// Wait for received interpolated values and store in GeoVaLs
  ASSERT(recvinterp_.size() == ntasks_);
  for (size_t jreq = 0; jreq < recv_req_.size(); ++jreq) {
    int ireq = -1;
    eckit::mpi::Status rst = comm_.waitAny(recv_req_, ireq);
    ASSERT(rst.error() == 0);
    size_t itask = recv_tasks_[ireq];
    ASSERT(itask >=0 && (size_t)itask < ntasks_);

    ASSERT(recvinterp_[itask].size() == myobs_index_by_task_[itask].size() * linsizes_);

    // Create non-owning views ("maps") into the interpolation results.
    const Eigen::Map<const Eigen::VectorX<size_t>> indices(myobs_index_by_task_[itask].data(),
                                                           myobs_index_by_task_[itask].size());
    // Each column contains the values of a single variable at a single level and all locations
    // with indices 'indices'. The columns are ordered first by level and then by variable.
    const Eigen::Map<const Eigen::MatrixXd> values(recvinterp_[itask].data(),
                                                   myobs_index_by_task_[itask].size(), linsizes_);

    size_t colOffset = 0;
    for (size_t jvar = 0; jvar < linvars_.size(); ++jvar) {
      const size_t numLevels = linvarsSizes_[jvar];
      geovals.fill(linvars_[jvar], indices, values.middleCols(colOffset, numLevels),
                   this->levelsTopDown_);
      colOffset += numLevels;
    }
    ASSERT(colOffset == linsizes_);
  }
  recv_req_.clear();
  recv_tasks_.clear();
  recvinterp_.clear();

// Clean-up send buffers (after making sure data has been sent)
  for (size_t jreq = 0; jreq < send_req_.size(); ++jreq) {
    int itask = -1;
    eckit::mpi::Status sst = comm_.waitAny(send_req_, itask);
    ASSERT(sst.error() == 0);
  }
  send_req_.clear();
  locinterp_.clear();
  geovalsTL_ = false;

  Log::trace() << "GetValues::fillGeoVaLsTL done" << std::endl;
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

  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
//  Mask obs outside time slot
    const util::TimeWindow timeSubWindow =
      timeWindow_.createSubWindow(dx.validTime(), hslot_);
    const std::vector<bool> mask =
      timeSubWindow.createTimeMask(obs_times_by_task_[jtask]);

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

  hslot_ = tstep/2;

// (Adjoint of) Send values interpolated locally (non-blocking)
// i.e. wait for receive of local sensitivities
  ASSERT(locinterp_.size() == ntasks_);
  for (size_t jreq = 0; jreq < send_req_.size(); ++jreq) {
    int itask = -1;
    eckit::mpi::Status sst = comm_.waitAny(send_req_, itask);
    ASSERT(sst.error() == 0);
    ASSERT(itask >=0 && (size_t)itask < ntasks_);
  }
  send_req_.clear();

// (Adjoint of) Allocate receive buffers and non blocking receive of interpolated values
// i.e. deallocate buffers (after making sure data has been sent)
  ASSERT(recvinterp_.size() == ntasks_);
  for (size_t jreq = 0; jreq < recv_req_.size(); ++jreq) {
    int itask = -1;
    eckit::mpi::Status rst = comm_.waitAny(recv_req_, itask);
    ASSERT(rst.error() == 0);
  }
  recv_req_.clear();
  recvinterp_.clear();

  Log::trace() << "GetValues::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::fillGeoVaLsAD(const GeoVaLs_ & geovals) {
  Log::trace() << "GetValues::fillGeoVaLsAD start" << std::endl;
  util::Timer timer("oops::GetValues", "fillGeoVaLsAD");

  const double missing = util::missingValue<double>();

// (Afjoint of) Clean-up send buffers
// i.e. allocate buffer and prepare to receive values
  ASSERT(locinterp_.empty());
  ASSERT(send_req_.empty());
  locinterp_.resize(ntasks_);
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    locinterp_[jtask].resize(obs_times_by_task_[jtask].size() * linsizes_, missing);
    if (locinterp_[jtask].size() > 0) {
      send_req_.push_back(comm_.iReceive(&locinterp_[jtask][0], locinterp_[jtask].size(),
                                         jtask, tag_));
    }
  }

// (Adjoint of) Wait for received interpolated values and store in GeoVaLs
// i.e. get values from GeoVaLs and send them
  ASSERT(recvinterp_.empty());
  ASSERT(recv_req_.empty());
  recvinterp_.resize(ntasks_);
  for (size_t jtask = 0; jtask < ntasks_; ++jtask) {
    const size_t nrecv = myobs_index_by_task_[jtask].size() * linsizes_;
    recvinterp_[jtask].resize(nrecv);

    if (nrecv > 0) {
      // Create non-owning views ("maps") into the interpolation results.
      const Eigen::Map<const Eigen::VectorX<size_t>> indices(myobs_index_by_task_[jtask].data(),
                                                             myobs_index_by_task_[jtask].size());
      // Each column contains the values of a single variable at a single level and all locations
      // with indices 'indices'. The columns are ordered first by level and then by variable.
      Eigen::Map<Eigen::MatrixXd> values(recvinterp_[jtask].data(),
                                         myobs_index_by_task_[jtask].size(), linsizes_);
      size_t colOffset = 0;
      for (size_t jvar = 0; jvar < linvars_.size(); ++jvar) {
        const size_t numLevels = linvarsSizes_[jvar];
        geovals.fillAD(linvars_[jvar], indices, values.middleCols(colOffset, numLevels),
                       this->levelsTopDown_);
        colOffset += numLevels;
      }
      ASSERT(colOffset == linsizes_);

      recv_req_.push_back(comm_.iSend(recvinterp_[jtask].data(), nrecv, jtask, tag_));
    }
  }

// Add MPI barrier to work around an intel MPI deadlock on some AMD platforms;
// barrier is added every N'th obs type
#ifdef INTELMPI_DEADLOCK_GETVALUES_LIMIT
  if ( tag_ % INTELMPI_DEADLOCK_GETVALUES_LIMIT == 0 ) comm_.barrier();
#endif

  Log::trace() << "GetValues::fillGeoVaLsAD" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
