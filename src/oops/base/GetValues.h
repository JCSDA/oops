/*
 * (C) Copyright 2020-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <Eigen/Core>

#include <algorithm>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Locations.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/LocalInterpolator.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Timer.h"
#include "oops/util/TimeWindow.h"

namespace oops {

/// \brief Fills GeoVaLs with requested variables at obs locations during model run
template <typename MODEL, typename OBS>
class GetValues : private util::ObjectCounter<GetValues<MODEL, OBS> > {
  typedef Geometry<MODEL>           Geometry_;
  typedef GeoVaLs<OBS>              GeoVaLs_;
  typedef Increment<MODEL>          Increment_;
  typedef LocalInterpolator<MODEL>  LocalInterpolator_;
  typedef Locations<OBS>            Locations_;
  typedef State<MODEL>              State_;

 public:
  static const std::string classname() {return "oops::GetValues";}

  GetValues(const eckit::Configuration &, const Geometry_ &,
            const util::TimeWindow &,
            const Locations_ &,
            const Variables &, const Variables & linvars = Variables());

  // Expose the LocalInterpolator's preprocess. This enables the user code to
  // call preprocess to process a State/Increment a single time even when multiple
  // GetValues are in use, thus optimizing certain algorithms.
  static void preprocess(State_ &);
  static void preprocess(Increment_ &);
  static void preprocessAD(Increment_ &);

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

/// Continuous DA update
  void updateGetVals(const eckit::Configuration &);

 private:
/// time-interpolation helper: adds contribution from this time to running total
  void incInterpValues(const util::DateTime &, const std::vector<bool> &,
                       int, int, const std::vector<double> &, std::vector<double>::iterator &);

  const eckit::mpi::Comm & comm_;
  const int ntasks_;
  int tag_;

  util::TimeWindow timeWindow_;
  util::Duration halfWidth_;
  bool requestLinearTimeInterpolation_ = false;  // request linear vs nearest-neighbor time interp
  bool doLinearTimeInterpolation_ = false;       // actually perform linear vs NN time interp

  const Variables geovars_;  // variables to interpolate and fill GeoVaLs with
  const Variables linvars_;

  // GetValues is responsible for populating the GeoVaLs associated with one
  // obs type (i.e., one ObsSpace), by interpolating model fields to observation
  // coordinates. Each obs type can be represented as a list of SamplingMethods,
  // where a sampling method pairs variables to coordinates. Using a list of
  // sampling methods is how one obs operator requests GeoVaLs that have different
  // variables at different coordinates.
  // To implement this, some GetValues data member range over the list of sampling
  // methods. "Sampling method" described by "sm" suffix.
  int nsms_ = 0;  // number of sampling methods
  std::vector<Variables> geovars_sm_;
  std::vector<Variables> linvars_sm_;
  std::vector<std::vector<size_t>> geolevels_sm_;
  std::vector<std::vector<size_t>> linlevels_sm_;
  std::vector<size_t> varsizes_sm_;
  std::vector<size_t> linsizes_sm_;

  // The overall interpolation matrix is broken up into sub-matrices, one per
  // each obs-owning task (i.e., task the interpolation result will be sent to)
  // per each sampling method. "Obs owning task" described by "ot" suffix.
  std::vector<std::vector<std::unique_ptr<LocalInterpolator_>>> interp_ot_sm_;
  // The obs times, grouped per each obs-owning task per each sampling method.
  std::vector<std::vector<std::vector<util::DateTime>>> times_ot_sm_;
  // The obs indices, grouped per each obs-interpolating task (i.e., the task
  // owning the model subdomain in which the observation lies, or equivalently
  // the task from which the interpolation result will be sent) per each
  // sampling method. "Model owning task" described by "mt" suffix.
  std::vector<std::vector<std::vector<size_t>>> indices_mt_sm_;

  std::vector<std::vector<double>> send_buffers_;
  std::vector<std::vector<double>> recv_buffers_;
  std::vector<eckit::mpi::Request> send_reqs_;
  std::vector<eckit::mpi::Request> recv_reqs_;
  std::vector<int> recv_tasks_;

  const bool levelsTopDown_;
  bool geovalsTL_ = false;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GetValues<MODEL, OBS>::GetValues(const eckit::Configuration & conf, const Geometry_ & geom,
                                 const util::TimeWindow & timeWindow,
                                 const Locations_ & locations,
                                 const Variables & vars, const Variables & linvars)
  : comm_(geom.getComm()), ntasks_(comm_.size()), tag_(789),
    timeWindow_(timeWindow), halfWidth_(0),
    geovars_(vars), linvars_(linvars),
    nsms_(locations.numSamplingMethods()),
    geovars_sm_(nsms_), linvars_sm_(nsms_),
    geolevels_sm_(nsms_), linlevels_sm_(nsms_),
    varsizes_sm_(nsms_), linsizes_sm_(nsms_),
    interp_ot_sm_(ntasks_), times_ot_sm_(ntasks_),
    indices_mt_sm_(ntasks_),
    send_buffers_(), recv_buffers_(), send_reqs_(), recv_reqs_(),
    levelsTopDown_(geom.levelsAreTopDown())
{
  Log::trace() << "GetValues::GetValues start" << std::endl;
  util::Timer timer("oops::GetValues", "GetValues");

  // each obs space creates a GetValues, each gets a unique MPI tag
  tag_ += this->created();

  // set the type of time-interpolation
  requestLinearTimeInterpolation_ = false;
  if (conf.has("time interpolation")) {
    const std::string value = conf.getString("time interpolation");
    if (value == "linear") {
      requestLinearTimeInterpolation_ = true;
    } else if (value != "nearest") {
      ABORT("GetValues::GetValues: time interpolation has an unsupported value.");
    }
  }

  for (auto var : geovars_) {
    geovars_sm_[locations.samplingMethodIndex(var)].push_back(var);
  }
  for (auto var : linvars_) {
    linvars_sm_[locations.samplingMethodIndex(var)].push_back(var);
  }

  for (int jsm = 0; jsm < nsms_; ++jsm) {
    geolevels_sm_[jsm] = geom.variableSizes(geovars_sm_[jsm]);
    linlevels_sm_[jsm] = geom.variableSizes(linvars_sm_[jsm]);
  }
  for (int jsm = 0; jsm < nsms_; ++jsm) {
    varsizes_sm_[jsm] = std::accumulate(
        geolevels_sm_[jsm].begin(), geolevels_sm_[jsm].end(), 0);
    linsizes_sm_[jsm] = std::accumulate(
        linlevels_sm_[jsm].begin(), linlevels_sm_[jsm].end(), 0);
  }

  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    interp_ot_sm_[jtask].resize(nsms_);
    times_ot_sm_[jtask].resize(nsms_);
    indices_mt_sm_[jtask].resize(nsms_);
  }

  const auto partitioner = LocalInterpolator_::makeTargetPartitioner(geom);

  std::vector<std::vector<double>> coords_mt(ntasks_);
  for (int jsm = 0; jsm < nsms_; ++jsm) {
    const auto & locs = locations.samplingMethod(jsm);
    const std::vector<double> & obslats = locs.latitudes();
    const std::vector<double> & obslons = locs.longitudes();
    const std::vector<util::DateTime> & obstimes = locs.times();

    // Assign obs to model-grid processors
    for (size_t jobs = 0; jobs < obstimes.size(); ++jobs) {
      const int itask = partitioner.interpolatingTask(obslats[jobs], obslons[jobs]);
      indices_mt_sm_[itask][jsm].push_back(jobs);
    }

    // Then allocate and copy
    for (int jtask = 0; jtask < ntasks_; ++jtask) {
      const size_t nb_obs_tsm = indices_mt_sm_[jtask][jsm].size();
      coords_mt[jtask].reserve(2 + 4 * nb_obs_tsm + coords_mt[jtask].size());
      coords_mt[jtask].push_back(jsm);
      coords_mt[jtask].push_back(nb_obs_tsm);
      for (size_t jobs = 0; jobs < nb_obs_tsm; ++jobs) {
        const size_t iobs = indices_mt_sm_[jtask][jsm][jobs];
        coords_mt[jtask].push_back(obslats[iobs]);
        coords_mt[jtask].push_back(obslons[iobs]);
        obstimes[iobs].serialize(coords_mt[jtask]);
      }
    }
  }

// Verify that an exception will not be thrown in the underlying eckit code.
  const std::string errorMessageSize =
    "The product of (the number of observation locations) and "
    "(the sum, over all model variables, of the number of levels in each GeoVaL) "
    "on this MPI rank is larger than the maximum integer, which will trigger an assertion "
    "in the underlying MPI communication code. One way to mitigate this is to reduce the "
    "number of locations on each MPI rank, either by increasing the number of "
    "ranks available or reducing the size of the input data set. "
    "Another way is to request fewer model variables.";
  const size_t maxSizeForCommunication = static_cast<size_t>(std::numeric_limits<int>::max());
  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    size_t biggestMessage = 0;
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      biggestMessage += indices_mt_sm_[jtask][jsm].size() * varsizes_sm_[jsm];
    }
    if (biggestMessage >= maxSizeForCommunication) {
      throw eckit::UserError(errorMessageSize, Here());
    }
  }

  // Exchange interpolation-target coordinates
  // from _mt (coords grouped by task doing the interp) to _ot (by obs owning task)
  std::vector<std::vector<double>> coords_ot(ntasks_);
  comm_.allToAll(coords_mt, coords_ot);

  // Setup interpolators
  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    size_t ii = 0;
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      const int expected_jsm = coords_ot[jtask][ii++];
      ASSERT(expected_jsm == jsm);
      const size_t nobs = coords_ot[jtask][ii++];
      std::vector<double> lats(nobs);
      std::vector<double> lons(nobs);
      times_ot_sm_[jtask][jsm].resize(nobs);
      for (size_t jobs = 0; jobs < nobs; ++jobs) {
        lats[jobs] = coords_ot[jtask][ii++];
        lons[jobs] = coords_ot[jtask][ii++];
        times_ot_sm_[jtask][jsm][jobs].deserialize(coords_ot[jtask], ii);
      }
      interp_ot_sm_[jtask][jsm] = std::make_unique<LocalInterpolator_>(conf, geom, lats, lons);
    }
    ASSERT(coords_ot[jtask].size() == ii);
  }

  Log::trace() << "GetValues::GetValues done" << std::endl;
}

// -----------------------------------------------------------------------------
//  Preprocess methods
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::preprocess(State_ & xx) {
  LocalInterpolator_::preprocess(xx);
}

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::preprocess(Increment_ & dx) {
  LocalInterpolator_::preprocess(dx);
}

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::preprocessAD(Increment_ & dx) {
  LocalInterpolator_::preprocessAD(dx);
}

// -----------------------------------------------------------------------------
//  Forward methods (called from nonlinear run)
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::initialize(const util::Duration & tstep) {
  Log::trace() << "GetValues::initialize start" << std::endl;
  util::Timer timer("oops::GetValues", "initialize");

  const double missing = util::missingValue<double>();

  ASSERT(send_buffers_.empty());
  send_buffers_.resize(ntasks_);
  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    size_t buffer_size = 0;
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      const size_t nobs = times_ot_sm_[jtask][jsm].size();
      buffer_size += nobs * varsizes_sm_[jsm];
    }
    send_buffers_[jtask].resize(buffer_size, missing);
  }

  // no need to do time interpolation if there is only one subwindow
  doLinearTimeInterpolation_ = requestLinearTimeInterpolation_ && (tstep < timeWindow_.length());
  halfWidth_ = doLinearTimeInterpolation_ ? tstep : tstep/2;

  Log::trace() << "GetValues::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::incInterpValues(
                    const util::DateTime & tCurrent, const std::vector<bool> & mask,
                    const int jtask,
                    const int jsm,
                    const std::vector<double> & tmp_buffer,
                    std::vector<double>::iterator & buffer_position)
{
  Log::trace() << "GetValues::incInterpValues start" << std::endl;

  const double missing = util::missingValue<double>();

// Get the state previous and next times and time-step
  const util::DateTime tPrevious = tCurrent - halfWidth_;
  const util::DateTime tNext = tCurrent + halfWidth_;
  const double dt = static_cast<double>(halfWidth_.toSeconds());

// Compute and add time weighted contribution from the input interpolated values
  auto & lhs_ptr = buffer_position;  // alias to a more helpful name within this scope
  auto rhs_ptr = tmp_buffer.begin();
  const int nb_obs = times_ot_sm_[jtask][jsm].size();
  for (size_t jf = 0; jf < geovars_sm_[jsm].size(); ++jf) {
    const int nb_levs = geolevels_sm_[jsm][jf];
    const auto shape = atlas::array::ArrayShape{nb_obs, nb_levs};

    // Get array views into interpolation results
    std::unique_ptr<atlas::array::Array> lhs_array(
        atlas::array::Array::wrap<double>(&*lhs_ptr, shape));
    auto lhs = atlas::array::make_view<double, 2>(*lhs_array);
    // Horrible cast to allow a standard ArrayView<double> from a const std::vector
    std::unique_ptr<atlas::array::Array> rhs_array(
        atlas::array::Array::wrap<double>(const_cast<double*>(&*rhs_ptr), shape));
    auto rhs = atlas::array::make_view<double, 2>(*rhs_array);

    for (int jloc = 0; jloc < nb_obs; ++jloc) {
      if (mask[jloc]) {
        // Compute time-interpolation weights
        const util::DateTime & obCurrentTime = times_ot_sm_[jtask][jsm][jloc];
        const bool isCurrentTime = (obCurrentTime == tCurrent);
        const bool isFirst = (obCurrentTime > tCurrent);
        double timeWeight = 0.;
        if (!isCurrentTime) {
          timeWeight = isFirst ?
            static_cast<double>((tNext - obCurrentTime).toSeconds())/dt :
            static_cast<double>((obCurrentTime - tPrevious).toSeconds())/dt;
        }

        for (int jlev = 0; jlev < nb_levs; ++jlev) {
          if (rhs(jloc, jlev) == missing) {
            lhs(jloc, jlev) = missing;
          } else if (isCurrentTime) {
            lhs(jloc, jlev) = rhs(jloc, jlev);
          } else if (isFirst) {
            lhs(jloc, jlev) = rhs(jloc, jlev) * timeWeight;
          } else if (lhs(jloc, jlev) != missing) {
            // Don't interpolate missing values
            lhs(jloc, jlev) += rhs(jloc, jlev) * timeWeight;
          }
        }
      }
    }

    const int step = nb_obs * nb_levs;
    std::advance(lhs_ptr, step);
    std::advance(rhs_ptr, step);
  }

  ASSERT(lhs_ptr == send_buffers_[jtask].end());
  ASSERT(rhs_ptr == tmp_buffer.end());

  Log::trace() << "GetValues::incInterpValues done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::process(const State_ & xx) {
  Log::trace() << "GetValues::process start" << std::endl;
  util::Timer timer("oops::GetValues", "process");

  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    const util::TimeWindow timeSubWindow =
      timeWindow_.createSubWindow(xx.validTime(), halfWidth_);

    auto buffer_position = send_buffers_[jtask].begin();

    for (int jsm = 0; jsm < nsms_; ++jsm) {
      // Mask obs outside time slot
      const std::vector<bool> mask =
        timeSubWindow.createTimeMask(times_ot_sm_[jtask][jsm]);

      const size_t nobs = times_ot_sm_[jtask][jsm].size();
      const size_t size = nobs * varsizes_sm_[jsm];

      // TODO(FH): refactor LocalInterpolator interface to write in-place into send_buffers_,
      //           using some sort of a view like gsl::span. This will avoid the allocation
      //           of tmp_buffer and copy from tmp_buffer to send_buffers_.
      std::vector<double> tmp_buffer(buffer_position,
                                     std::next(buffer_position, size));
      interp_ot_sm_[jtask][jsm]->apply(geovars_sm_[jsm], xx, mask, tmp_buffer);

      if (doLinearTimeInterpolation_) {
        // call below will advance buffer_position
        incInterpValues(xx.validTime(), mask, jtask, jsm, tmp_buffer, buffer_position);
      } else {
        std::copy(tmp_buffer.begin(), tmp_buffer.end(), buffer_position);
        std::advance(buffer_position, size);
      }
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
  ASSERT(recv_buffers_.empty());
  ASSERT(recv_reqs_.empty());
  ASSERT(recv_tasks_.empty());
  recv_buffers_.resize(ntasks_);
  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    size_t buffer_size = 0;
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      const size_t nobs = indices_mt_sm_[jtask][jsm].size();
      buffer_size += nobs * varsizes_sm_[jsm];
    }
    if (buffer_size > 0) {
      recv_buffers_[jtask].resize(buffer_size);
      recv_reqs_.push_back(comm_.iReceive(&recv_buffers_[jtask][0], buffer_size, jtask, tag_));
      recv_tasks_.push_back(jtask);
    }
  }

// Send values interpolated locally (non-blocking)
  ASSERT(send_reqs_.empty());
  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    if (send_buffers_[jtask].size() > 0) {
      send_reqs_.push_back(comm_.iSend(&send_buffers_[jtask][0], send_buffers_[jtask].size(),
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
  ASSERT(recv_buffers_.size() == static_cast<size_t>(ntasks_));
  for (size_t jreq = 0; jreq < recv_reqs_.size(); ++jreq) {
    int ireq = -1;
    eckit::mpi::Status rst = comm_.waitAny(recv_reqs_, ireq);
    ASSERT(rst.error() == 0);
    const int itask = recv_tasks_[ireq];
    ASSERT(itask >=0 && itask < ntasks_);

    size_t expected_size = 0;
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      expected_size += indices_mt_sm_[itask][jsm].size() * varsizes_sm_[jsm];
    }
    ASSERT(recv_buffers_[itask].size() == expected_size);

    auto offset = recv_buffers_[itask].begin();
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      // Create non-owning views ("maps") into the interpolation results.
      const size_t numLocs = indices_mt_sm_[itask][jsm].size();
      const Eigen::Map<const Eigen::VectorX<size_t>> indices(
          indices_mt_sm_[itask][jsm].data(), numLocs);

      for (size_t jvar = 0; jvar < geovars_sm_[jsm].size(); ++jvar) {
        const size_t numLevels = geolevels_sm_[jsm][jvar];
        // View the buffer as an Eigen matrix; the contiguous dimension ranges over
        // the height of the columns, the strided dimension ranges over the obs
        // locations with indices `indices`.
        const Eigen::Map<const Eigen::MatrixXd> values(&*offset, numLevels, numLocs);
        geovals.fill(geovars_sm_[jsm][jvar], indices, values, levelsTopDown_);
        std::advance(offset, numLevels * numLocs);
      }
    }
    ASSERT(offset == recv_buffers_[itask].end());
  }
  recv_reqs_.clear();
  recv_tasks_.clear();
  recv_buffers_.clear();

// Clean-up send buffers (after making sure data has been sent)
  for (size_t jreq = 0; jreq < send_reqs_.size(); ++jreq) {
    int itask = -1;
    eckit::mpi::Status sst = comm_.waitAny(send_reqs_, itask);
    ASSERT(sst.error() == 0);
  }
  send_reqs_.clear();
  send_buffers_.clear();

  Log::trace() << "GetValues::fillGeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------
//  TL methods
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::initializeTL(const util::Duration & tstep) {
  Log::trace() << "GetValues::initializeTL start" << std::endl;
  util::Timer timer("oops::GetValues", "initializeTL");

  const double missing = util::missingValue<double>();

  ASSERT(send_buffers_.empty());
  send_buffers_.resize(ntasks_);
  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    size_t buffer_size = 0;
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      const size_t nobs = times_ot_sm_[jtask][jsm].size();
      buffer_size += nobs * linsizes_sm_[jsm];
    }
    send_buffers_[jtask].resize(buffer_size, missing);
  }
  halfWidth_ = tstep/2;
  geovalsTL_ = true;

  Log::trace() << "GetValues::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::processTL(const Increment_ & dx) {
  Log::trace() << "GetValues::processTL start" << std::endl;
  util::Timer timer("oops::GetValues", "processTL");

  for (int jtask = 0; jtask < ntasks_; ++jtask) {
//  Mask obs outside time slot
    const util::TimeWindow timeSubWindow =
      timeWindow_.createSubWindow(dx.validTime(), halfWidth_);

    auto buffer_position = send_buffers_[jtask].begin();

    for (int jsm = 0; jsm < nsms_; ++jsm) {
      const std::vector<bool> mask =
        timeSubWindow.createTimeMask(times_ot_sm_[jtask][jsm]);

      const size_t nobs = times_ot_sm_[jtask][jsm].size();
      const size_t size = nobs * linsizes_sm_[jsm];

      // TODO(FH): refactor LocalInterpolator interface to write in-place into send_buffers_,
      //           using some sort of a view like gsl::span. This will avoid the allocation
      //           of tmp_buffer and copy from tmp_buffer to send_buffers_.
      std::vector<double> tmp_buffer(buffer_position,
                                     std::next(buffer_position, size));
      interp_ot_sm_[jtask][jsm]->apply(linvars_sm_[jsm], dx, mask, tmp_buffer);
      std::copy(tmp_buffer.begin(), tmp_buffer.end(), buffer_position);
      std::advance(buffer_position, size);
    }
  }

  Log::trace() << "GetValues::processTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::finalizeTL() {
  Log::trace() << "GetValues::finalizeTL start" << std::endl;
  util::Timer timer("oops::GetValues", "finalizeTL");

// Allocate receive buffers and non blocking receive of interpolated values
  ASSERT(recv_buffers_.empty());
  ASSERT(recv_reqs_.empty());
  ASSERT(recv_tasks_.empty());
  recv_buffers_.resize(ntasks_);
  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    size_t buffer_size = 0;
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      const size_t nobs = indices_mt_sm_[jtask][jsm].size();
      buffer_size += nobs * linsizes_sm_[jsm];
    }
    if (buffer_size > 0) {
      recv_buffers_[jtask].resize(buffer_size);
      recv_reqs_.push_back(comm_.iReceive(&recv_buffers_[jtask][0], buffer_size, jtask, tag_));
      recv_tasks_.push_back(jtask);
    }
  }

// Send values interpolated locally (non-blocking)
  ASSERT(send_reqs_.empty());
  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    if (send_buffers_[jtask].size() > 0) {
      send_reqs_.push_back(comm_.iSend(&send_buffers_[jtask][0], send_buffers_[jtask].size(),
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
  ASSERT(recv_buffers_.size() == static_cast<size_t>(ntasks_));
  for (size_t jreq = 0; jreq < recv_reqs_.size(); ++jreq) {
    int ireq = -1;
    eckit::mpi::Status rst = comm_.waitAny(recv_reqs_, ireq);
    ASSERT(rst.error() == 0);
    const int itask = recv_tasks_[ireq];
    ASSERT(itask >=0 && itask < ntasks_);

    size_t expected_size = 0;
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      expected_size += indices_mt_sm_[itask][jsm].size() * linsizes_sm_[jsm];
    }
    ASSERT(recv_buffers_[itask].size() == expected_size);

    auto offset = recv_buffers_[itask].begin();
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      // Create non-owning views ("maps") into the interpolation results.
      const size_t numLocs = indices_mt_sm_[itask][jsm].size();
      const Eigen::Map<const Eigen::VectorX<size_t>> indices(
          indices_mt_sm_[itask][jsm].data(), numLocs);

      for (size_t jvar = 0; jvar < linvars_sm_[jsm].size(); ++jvar) {
        const size_t numLevels = linlevels_sm_[jsm][jvar];
        // View the buffer as an Eigen matrix; the contiguous dimension ranges over
        // the height of the columns, the strided dimension ranges over the obs
        // locations with indices `indices`.
        const Eigen::Map<const Eigen::MatrixXd> values(&*offset, numLevels, numLocs);
        geovals.fill(linvars_sm_[jsm][jvar], indices, values, levelsTopDown_);
        std::advance(offset, numLevels * numLocs);
      }
    }
    ASSERT(offset == recv_buffers_[itask].end());
  }
  recv_reqs_.clear();
  recv_tasks_.clear();
  recv_buffers_.clear();

// Clean-up send buffers (after making sure data has been sent)
  for (size_t jreq = 0; jreq < send_reqs_.size(); ++jreq) {
    int itask = -1;
    eckit::mpi::Status sst = comm_.waitAny(send_reqs_, itask);
    ASSERT(sst.error() == 0);
  }
  send_reqs_.clear();
  send_buffers_.clear();
  geovalsTL_ = false;

  Log::trace() << "GetValues::fillGeoVaLsTL done" << std::endl;
}

// -----------------------------------------------------------------------------
//  AD methods
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::initializeAD() {
  Log::trace() << "GetValues::initializeAD start" << std::endl;
  util::Timer timer("oops::GetValues", "initializeAD");
  send_buffers_.clear();
  Log::trace() << "GetValues::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::processAD(Increment_ & dx) {
  Log::trace() << "GetValues::processAD start" << std::endl;
  util::Timer timer("oops::GetValues", "processAD");

  for (int jtask = 0; jtask < ntasks_; ++jtask) {
//  Mask obs outside time slot
    const util::TimeWindow timeSubWindow =
      timeWindow_.createSubWindow(dx.validTime(), halfWidth_);

    auto buffer_position = send_buffers_[jtask].begin();

    for (int jsm = 0; jsm < nsms_; ++jsm) {
      const std::vector<bool> mask =
        timeSubWindow.createTimeMask(times_ot_sm_[jtask][jsm]);

      const size_t nobs = times_ot_sm_[jtask][jsm].size();
      const size_t size = nobs * linsizes_sm_[jsm];

      // TODO(FH): refactor LocalInterpolator interface to write in-place into send_buffers_,
      //           using some sort of a view like gsl::span. This will avoid the allocation
      //           of tmp_buffer and copy from send_buffers_ to tmp_buffer.
      const std::vector<double> tmp_buffer(buffer_position,
                                           std::next(buffer_position, size));
      interp_ot_sm_[jtask][jsm]->applyAD(linvars_sm_[jsm], dx, mask, tmp_buffer);
      std::advance(buffer_position, size);
    }
  }

  Log::trace() << "GetValues::processAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::finalizeAD(const util::Duration & tstep) {
  Log::trace() << "GetValues::finalizeAD start" << std::endl;
  util::Timer timer("oops::GetValues", "finalizeAD");

  halfWidth_ = tstep/2;

// (Adjoint of) Send values interpolated locally (non-blocking)
// i.e. wait for receive of local sensitivities
  ASSERT(send_buffers_.size() == static_cast<size_t>(ntasks_));
  for (size_t jreq = 0; jreq < send_reqs_.size(); ++jreq) {
    int itask = -1;
    eckit::mpi::Status sst = comm_.waitAny(send_reqs_, itask);
    ASSERT(sst.error() == 0);
    ASSERT(itask >=0 && itask < ntasks_);
  }
  send_reqs_.clear();

// (Adjoint of) Allocate receive buffers and non blocking receive of interpolated values
// i.e. deallocate buffers (after making sure data has been sent)
  ASSERT(recv_buffers_.size() == static_cast<size_t>(ntasks_));
  for (size_t jreq = 0; jreq < recv_reqs_.size(); ++jreq) {
    int itask = -1;
    eckit::mpi::Status rst = comm_.waitAny(recv_reqs_, itask);
    ASSERT(rst.error() == 0);
  }
  recv_reqs_.clear();
  recv_buffers_.clear();

  Log::trace() << "GetValues::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::fillGeoVaLsAD(const GeoVaLs_ & geovals) {
  Log::trace() << "GetValues::fillGeoVaLsAD start" << std::endl;
  util::Timer timer("oops::GetValues", "fillGeoVaLsAD");

  const double missing = util::missingValue<double>();

// (Adjoint of) Clean-up send buffers
// i.e. allocate buffer and prepare to receive values
  ASSERT(send_buffers_.empty());
  ASSERT(send_reqs_.empty());
  send_buffers_.resize(ntasks_);
  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    size_t buffer_size = 0;
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      const size_t nobs = times_ot_sm_[jtask][jsm].size();
      buffer_size += nobs * linsizes_sm_[jsm];
    }
    if (buffer_size > 0) {
      send_buffers_[jtask].resize(buffer_size, missing);
      send_reqs_.push_back(comm_.iReceive(&send_buffers_[jtask][0], buffer_size, jtask, tag_));
    }
  }

// (Adjoint of) Wait for received interpolated values and store in GeoVaLs
// i.e. get values from GeoVaLs and send them
  ASSERT(recv_buffers_.empty());
  ASSERT(recv_reqs_.empty());
  recv_buffers_.resize(ntasks_);
  for (int jtask = 0; jtask < ntasks_; ++jtask) {
    size_t buffer_size = 0;
    for (int jsm = 0; jsm < nsms_; ++jsm) {
      const size_t nobs = indices_mt_sm_[jtask][jsm].size();
      buffer_size += nobs * linsizes_sm_[jsm];
    }
    if (buffer_size > 0) {
      recv_buffers_[jtask].resize(buffer_size);

      auto offset = recv_buffers_[jtask].begin();
      for (int jsm = 0; jsm < nsms_; ++jsm) {
        // Create non-owning views ("maps") into the interpolation results.
        const size_t numLocs = indices_mt_sm_[jtask][jsm].size();
        const Eigen::Map<const Eigen::VectorX<size_t>> indices(
            indices_mt_sm_[jtask][jsm].data(), numLocs);

        for (size_t jvar = 0; jvar < linvars_sm_[jsm].size(); ++jvar) {
          const size_t numLevels = linlevels_sm_[jsm][jvar];
          // View the buffer as an Eigen matrix; the contiguous dimension ranges over
          // the height of the columns, the strided dimension ranges over the obs
          // locations with indices `indices`.
          Eigen::Map<Eigen::MatrixXd> values(&*offset, numLevels, numLocs);
          geovals.fillAD(linvars_sm_[jsm][jvar], indices, values, levelsTopDown_);
          std::advance(offset, numLevels * numLocs);
        }
      }
      ASSERT(offset == recv_buffers_[jtask].end());

      recv_reqs_.push_back(comm_.iSend(recv_buffers_[jtask].data(), buffer_size, jtask, tag_));
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

template <typename MODEL, typename OBS>
void GetValues<MODEL, OBS>::updateGetVals(const eckit::Configuration & cdaConfig) {
  if (cdaConfig.has("time window")) {
    timeWindow_ = util::TimeWindow(cdaConfig.getSubConfiguration("time window"));
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
