/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_LOCALIZATION_H_
#define OOPS_BASE_LOCALIZATION_H_

#include <Eigen/Dense>

#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/exception/Exceptions.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/generic/LocalizationBase.h"
#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// \brief Abstract model-space localization class used by high level algorithms
/// and applications.
///
/// Note: to see methods that need to be implemented in a generic Localization
/// implementation, see LocalizationBase class in generic/LocalizationBase.h.
/// To see methods that need to be implemented in a MODEL-specific Localization
/// implementation, see interface::LocalizationBase class in
/// interface/LocalizationBase.h.
template <typename MODEL>
class Localization : public util::Printable,
                     private boost::noncopyable,
                     private util::ObjectCounter<Localization<MODEL> > {
  typedef Geometry<MODEL>     Geometry_;
  typedef Increment<MODEL>    Increment_;
  typedef Increment4D<MODEL>  Increment4D_;
  typedef LocalizationBase<MODEL>   LocBase_;

 public:
  static const std::string classname() {return "oops::Localization";}

  /// Set up Localization for \p geometry, configured with \p conf
  Localization(const Geometry_ & geometry,
               const oops::Variables & incVars,
               const eckit::Configuration & conf);
  ~Localization();

  /// Randomize \p dx and apply 4D localization. All 3D blocks of the 4D localization
  /// matrix are the same (and defined by 3D localization loc_)
  virtual void randomize(Increment4D_ & dx) const;
  /// Apply 4D localization. All 3D blocks of the 4D localization matrix are the same
  /// (and defined by 3D localization loc_)
  virtual void multiply(Increment4D_ & dx) const;

 private:
  /// Print, used in logging
  void print(std::ostream &) const override;

  std::unique_ptr<util::Timer> timeConstr_;

  /// Pointer to the Localization implementation
  std::unique_ptr<LocBase_> loc_;

  // Communication method (standard, fast, aggressive)
  // standard: reproducible, axpy based on increment interfaces
  // fast: reproducible, axpy performed on serialized vectors
  // aggressive (only in the case without time decay):
  //     non reproducible (to machine precision), sums by mpi reduction
  const std::string commMode_;

  // Time communicator and rank
  const eckit::mpi::Comm & comm_;
  const size_t ntimes_;
  const size_t mytime_;

  // Time decay switch
  bool hasTimeDecay_;

  // Lower factor of the time decay matrix
  Eigen::MatrixXd TDLower_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::Localization(const Geometry_ & geometry,
                                  const oops::Variables & incVars,
                                  const eckit::Configuration & conf)
  : timeConstr_(new util::Timer(classname(), "Localization")),
    loc_(),
    commMode_(conf.getString("communication mode", "standard")),
    comm_(geometry.timeComm()),
    ntimes_(geometry.timeComm().size()),
    mytime_(geometry.timeComm().rank())
{
  Log::trace() << "Localization<MODEL>::Localization starting" << std::endl;

  // Define whether localization should be initialized on this task
  bool initLocOnMyTask(false);

  // Initialize time-decay flag
  hasTimeDecay_ = false;

  if (conf.has("time decay")) {
    // Duplicated localization with time decay
    initLocOnMyTask = true;

    // Set time-decay flag
    hasTimeDecay_ = true;

    // Check communication mode
    if (!(commMode_ == "standard" || commMode_ == "fast")) {
      ABORT("wrong communication option for this localization: " + commMode_);
    }

    // Get dates of each sub-window
    std::vector<util::DateTime> dates = {util::DateTime(conf.getString("date"))};
    oops::mpi::allGatherv(comm_, dates);

    // Get time-decay duration
    const double timeDecay = static_cast<double>(util::Duration(conf.getString("time decay"))
      .toSeconds());

    // Time decay matrix
    Eigen::MatrixXd TD(ntimes_, ntimes_);
    for (size_t i=0; i < ntimes_; ++i) {
      for (size_t j=0; j < ntimes_; ++j) {
        if (timeDecay > 0.0) {
          // Duration between subwindows
          const double diff = static_cast<double>((dates[i]-dates[j]).toSeconds());

          // Gaussian decay
          TD(i, j) = std::exp(-0.5*(diff*diff)/(timeDecay*timeDecay));
        } else {
          // Cut-off cross-localization
          if (i == j) {
            TD(i, j) = 1.0;
          } else {
            TD(i, j) = 0.0;
          }
        }
      }
    }

    // Cholesky decomposition to get the time factor
    Eigen::LDLT<Eigen::MatrixXd> ldlt(TD);
    Eigen::MatrixXd L(ntimes_, ntimes_);
    L = ldlt.matrixL();
    Eigen::ArrayXXd d(ntimes_, 1);
    d =  ldlt.vectorD().array();
    for (size_t i = 0; i < ntimes_; i++) {
      if (d(i, 0) < 0.0) d(i, 0) = 0.0;
    }
    Eigen::MatrixXd sqrtd(ntimes_, 1);
    sqrtd = d.sqrt();
    TDLower_ = ldlt.transpositionsP().transpose()*L*sqrtd.asDiagonal();
    Eigen::MatrixXd TDerror(ntimes_, ntimes_);
    TDerror = (TDLower_*(TDLower_.transpose())-TD).cwiseAbs();
    Log::info() << "Time-decay matrix: " << std::endl << TD << std::endl;
    Log::info() << "LDLT matrix error: " << TDerror.maxCoeff() << std::endl;
  } else {
    // Check communication mode
    if (!(commMode_ == "standard" || commMode_ == "fast" || commMode_ == "aggressive")) {
      ABORT("wrong communication option for this localization: " + commMode_);
    }

    // Duplicated localization without time decay
    initLocOnMyTask = (mytime_ == 0);
  }

  // Initialize localization
  if (initLocOnMyTask) {
    loc_ = std::move(LocalizationFactory<MODEL>::create(geometry, incVars, conf));
  }

  Log::trace() << "Localization<MODEL>::Localization done" << std::endl;
  timeConstr_.reset();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::~Localization() {
  Log::trace() << "Localization<MODEL>::~Localization starting" << std::endl;
  util::Timer timer(classname(), "~Localization");
  loc_.reset();
  Log::trace() << "Localization<MODEL>::~Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Localization<MODEL>::randomize(Increment4D_ & dx) const {
  Log::trace() << "Localization<MODEL>::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");

  if (dx.local_time_size() > 1)
    throw eckit::NotImplemented("4D localization not fully implemented", Here());

  // Save original time value
  util::DateTime tsub = dx[0].validTime();

  if (hasTimeDecay_) {
    // Duplicated localization with time-decay

    // Set time to 0
    util::DateTime t0(0, 0);
    dx[0].updateTime(t0-dx[0].validTime());

    // Apply 3D localization
    loc_->randomize(dx[0]);

    // Create temporary increment
    Increment_ dxtmp(dx[0], true);

    // Set output increment to zero
    dx.zero();

    // Apply lower triangular matrix of weights
    if (commMode_ == "standard") {
      // Apply lower triangular matrix of weights
      for (size_t j=0; j <= mytime_; ++j) {
        dx[0].axpy(TDLower_(mytime_, mytime_-j), dxtmp, false);
        size_t dest = mytime_ + 1;
        if (mytime_ == ntimes_ - 1 ) dest = comm_.procNull();
        size_t src = mytime_ - 1;
        if (mytime_ == j) src = comm_.procNull();
        oops::mpi::sendReceiveReplace(comm_, dxtmp, dest, 0, src, 0);
      }
    } else if (commMode_ == "fast") {
      // Serialize the output
      std::vector<double> dx_s;
      dx[0].serialize(dx_s);
      size_t sz = dx[0].serialSize();
      Eigen::Map<Eigen::VectorXd> dx_v(dx_s.data(), sz);

      // Serialize the temporary increment
      std::vector<double> dxtmp_s;
      dxtmp.serialize(dxtmp_s);
      Eigen::Map<Eigen::VectorXd> dxtmp_v(dxtmp_s.data(), sz);

      // Apply lower triangular matrix of weights on serialized vector
      for (size_t j=0; j <= mytime_; ++j) {
        dx_v += TDLower_(mytime_, mytime_-j) * dxtmp_v;
        size_t dest = mytime_ + 1;
        if (mytime_ == ntimes_ - 1 ) dest = comm_.procNull();
        size_t src = mytime_ - 1;
        if (mytime_ == j) src = comm_.procNull();
        eckit::mpi::Status status = comm_.sendReceiveReplace(dxtmp_s.data(), sz,
                                                             dest, 0, src, 0);
        }

      // Deserialize and store the result
      size_t ii = 0;
      dx[0].deserialize(dx_s, ii);
    }
  } else {
    // Duplicated localization without time decay
    if (mytime_ == 0) {
      // Apply 3D localization
      loc_->randomize(dx[0]);
    }

    // Broadcast
    oops::mpi::broadcast(comm_, dx[0], 0);
  }

  // Set time back to original value
  dx[0].updateTime(tsub - dx[0].validTime());

  Log::trace() << "Localization<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Localization<MODEL>::multiply(Increment4D_ & dx) const {
  Log::trace() << "Localization<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");

  // Save original time value
  util::DateTime tsub = dx[0].validTime();

  if (hasTimeDecay_) {
    if (dx.local_time_size() > 1)
      throw eckit::NotImplemented("4D localization not fully implemented", Here());
    // Duplicated localization with time-decay

    // Set time to 0
    util::DateTime t0(0, 0);
    dx[0].updateTime(t0-dx[0].validTime());

    // Create temporary increment
    Increment_ dxtmp(dx[0]);

    if (commMode_ == "fast") {
      // Serialize the input
      std::vector<double> dx_s;
      dx[0].serialize(dx_s);
      size_t sz = dx[0].serialSize();
      Eigen::Map<Eigen::VectorXd> dx_v(dx_s.data(), sz);

      // Set temporary increment to zero and serialize it
      dxtmp.zero();
      std::vector<double> dxtmp_s;
      dxtmp.serialize(dxtmp_s);
      Eigen::Map<Eigen::VectorXd> dxtmp_v(dxtmp_s.data(), sz);

      // Apply upper triangular matrix of weights on serialized vectors
      for (size_t j=0; j < ntimes_-mytime_; ++j) {
        dxtmp_v += TDLower_(mytime_+j, mytime_) * dx_v;
        size_t dest = mytime_ - 1;
        if (mytime_ == 0) dest = comm_.procNull();
        size_t src = mytime_ + 1;
        if (mytime_ == ntimes_ - 1 - j) src = comm_.procNull();
        eckit::mpi::Status status = comm_.sendReceiveReplace(dx_s.data(), sz,
                                                             dest, 0, src, 0);
      }

      // Apply 3D localization
      size_t ii = 0;
      dxtmp.deserialize(dxtmp_s, ii);
      dxtmp.updateTime(tsub - dxtmp.validTime());

      loc_->multiply(dxtmp);

      dxtmp.updateTime(t0-dxtmp.validTime());
      dxtmp_s.clear();
      dxtmp.serialize(dxtmp_s);

      // Set output increment to zero
      dx[0].zero();
      dx_s.clear();
      dx[0].serialize(dx_s);

      // Apply lower triangular matrix of weights on serialized vector
      for (size_t j=0; j <= mytime_; ++j) {
        dx_v += TDLower_(mytime_, mytime_-j) * dxtmp_v;
        size_t dest = mytime_ + 1;
        if (mytime_ == ntimes_ - 1) dest = comm_.procNull();
        size_t src = mytime_ - 1;
        if (mytime_ == j) src = comm_.procNull();
        eckit::mpi::Status status = comm_.sendReceiveReplace(dxtmp_s.data(), sz,
                                                             dest, 0, src, 0);
      }

      // Deserialize and store the result
      ii = 0;
      dx[0].deserialize(dx_s, ii);
    } else if (commMode_ == "standard") {
      // Set temporary increment to zero
      dxtmp.zero();

      // Apply upper triangular matrix of weights
      for (size_t j=0; j < ntimes_-mytime_; ++j) {
        dxtmp.axpy(TDLower_(mytime_+j, mytime_), dx[0], false);
        size_t dest = mytime_ - 1;
        if (mytime_ == 0) dest = comm_.procNull();
        size_t src = mytime_ + 1;
        if (mytime_ == ntimes_ - 1 - j) src = comm_.procNull();
        oops::mpi::sendReceiveReplace(comm_, dx[0], dest, 0, src, 0);
      }

      // Apply 3D localization
      loc_->multiply(dxtmp);

      // Set output increment to zero
      dx[0].zero();

      // Apply lower triangular matrix of weights
      for (size_t j=0; j <= mytime_; ++j) {
        dx[0].axpy(TDLower_(mytime_, mytime_-j), dxtmp, false);
        size_t dest = mytime_ + 1;
        if (mytime_ == ntimes_ - 1) dest = comm_.procNull();
        size_t src = mytime_ - 1;
        if (mytime_ == j) src = comm_.procNull();
        oops::mpi::sendReceiveReplace(comm_, dxtmp, dest, 0, src, 0);
      }
    }

    // Set time back to original value
    dx[0].updateTime(tsub - dx[0].validTime());

  } else {
  // Use Mark Buehner's trick to save CPU when applying the same 3D localization for all
  // 3D blocks of the 4D localization matrix:
  // L_4D = ( L_3D L_3D L_3D ) = ( Id ) L_3D ( Id Id Id )
  //        ( L_3D L_3D L_3D )   ( Id )
  //        ( L_3D L_3D L_3D )   ( Id )
  // so if :
  // x_4D = ( x_1 )
  //        ( x_2 )
  //        ( x_3 )
  // then:
  // L_4D x_4D = (Id) L_3D (Id Id Id) (x_1) = (Id) L_3D (x_1+x_2+x_3) = (L_3D ( x_1 + x_2 + x_3 ))
  //             (Id)                 (x_2)   (Id)                      (L_3D ( x_1 + x_2 + x_3 ))
  //             (Id)                 (x_3)   (Id)                      (L_3D ( x_1 + x_2 + x_3 ))
  // Reference in section 3.4.2. of https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.2325.

    if (commMode_ == "standard") {
      std::vector<util::DateTime> times(dx.validTimes());
      for (size_t jt = 1; jt < dx.local_time_size(); ++jt) {
        dx[0].axpy(1.0, dx[jt], false);
      }

      if (dx.commTime().rank() > 0) {
        oops::mpi::send(dx.commTime(), dx[0], 0, 0);
      } else {
        Increment_ dxtmp(dx[0], false);
        for (size_t jj = 1; jj < dx.commTime().size(); ++jj) {
          oops::mpi::receive(comm_, dxtmp, jj, 0);
          dx[0].axpy(1.0, dxtmp, false);
        }

        // Apply 3D localization
        loc_->multiply(dx[0]);
      }
      // Broadcast
      oops::mpi::broadcast(dx.commTime(), dx[0], 0);

      dx[0].updateTime(times[0] - dx[0].validTime());
      for (size_t jt = 1; jt < dx.local_time_size(); ++jt) {
        dx[jt] = dx[0];
        dx[jt].updateTime(times[jt] - dx[jt].validTime());
      }

    } else if (commMode_ == "fast") {
      if (dx.local_time_size() > 1)
        throw eckit::NotImplemented("4D localization not fully implemented", Here());
      // Set time to 0
      util::DateTime t0(0, 0);
      dx[0].updateTime(t0-dx[0].validTime());

      // Serialize the input
      std::vector<double> dx_s;
      dx[0].serialize(dx_s);
      size_t sz = dx[0].serialSize();
      Eigen::Map<Eigen::VectorXd> dx_v(dx_s.data(), sz);

      if (mytime_ > 0) {
        // Send to root task
        comm_.send(dx_s.data(), sz, 0, 0);
      } else {
        // Sum over timeslots
        Eigen::VectorXd dxtmp_v(sz);
        eckit::mpi::Status status;
        for (size_t jj = 1; jj < ntimes_; ++jj) {
          status = comm_.receive(dxtmp_v.data(), sz, static_cast<int>(jj), 0);
          dx_v += dxtmp_v;
        }

        // Apply 3D localization
        size_t ii = 0;
        dx[0].deserialize(dx_s, ii);
        dx[0].updateTime(tsub - dx[0].validTime());
        loc_->multiply(dx[0]);
      }
      // Broadcast
      oops::mpi::broadcast(comm_, dx[0], 0);

      if (mytime_ > 0) {
        // Set time back to original value
        dx[0].updateTime(tsub - dx[0].validTime());
      }

    } else if (commMode_ == "aggressive") {
      if (dx.local_time_size() > 1)
        throw eckit::NotImplemented("4D localization not fully implemented", Here());
      if (mytime_ > 0) {
        // Set time to 0
        util::DateTime t0(0, 0);
        dx[0].updateTime(t0-dx[0].validTime());
      }

      // Reduce on mytime_ 0
      oops::mpi::reduceInPlace(comm_, dx[0], 0);

      if (mytime_ == 0) {
        // Apply 3D localization
        loc_->multiply(dx[0]);
      }

      // Broadcast
      oops::mpi::broadcast(comm_, dx[0], 0);

      if (mytime_ > 0) {
        // Set time back to original value
        dx[0].updateTime(tsub - dx[0].validTime());
      }
    }
  }

  Log::trace() << "Localization<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Localization<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Localization<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *loc_;
  Log::trace() << "Localization<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_LOCALIZATION_H_
