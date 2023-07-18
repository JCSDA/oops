/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/gatherPrint.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template<typename DATA> class DataSetBase : public util::Printable {
 public:
  virtual ~DataSetBase() = default;

  DataSetBase & operator=(const DataSetBase &);

  const Variables & variables() const {this->check_consistency(); return dataset_[0].variables();}
  const std::vector<util::DateTime> validTimes() const;

  bool is_3d() const {return nmembers_ == 1 && ntimes_ == 1;}
  bool is_4d() const {return nmembers_ == 1;}
  size_t ens_size() const {return nmembers_;}
  size_t time_size() const {return ntimes_;}
  const std::vector<util::DateTime> & times() const {return times_;}

  DATA & operator()(const size_t it, const size_t im) {return this->data(it, im);}
  const DATA & operator()(const size_t it, const size_t im) const {return this->data(it, im);}

  size_t size() const {return dataset_.size();}
  DATA & operator[](const int ii) {return dataset_[ii];}
  const DATA & operator[](const int ii) const {return dataset_[ii];}

  void shift_forward();
  void shift_backward();

 protected:
  DataSetBase(const std::vector<util::DateTime> &, const eckit::mpi::Comm &,
              const std::vector<int> &, const eckit::mpi::Comm &);
  DataSetBase(const eckit::mpi::Comm &, const eckit::mpi::Comm &);
  DataSetBase(const DataSetBase &) = default;

  void check_consistency() const;
  void check_consistency(const DataSetBase &) const;

  DATA & data(const size_t it, const size_t im) {return dataset_.at(im * localtimes_ + it);}
  const DATA & data(const size_t it, const size_t im) const
                                                {return dataset_.at(im * localtimes_ + it);}

  size_t ntimes_;
  size_t localtimes_;
  std::vector<util::DateTime> times_;
  const eckit::mpi::Comm & commTime_;
  util::Duration subWinLength_;

  size_t nmembers_;
  size_t localmembers_;
  std::vector<int> members_;
  const eckit::mpi::Comm & commEns_;
  std::vector<size_t> mymembers_;

  std::vector<DATA> dataset_;

 private:
  void print(std::ostream &) const;
  virtual std::string classname() const = 0;
};

// -----------------------------------------------------------------------------

template<typename DATA>
DataSetBase<DATA>::DataSetBase(const std::vector<util::DateTime> & times,
                               const eckit::mpi::Comm & commTime,
                               const std::vector<int> & members,
                               const eckit::mpi::Comm & commEns)
  : ntimes_(times.size()), localtimes_(), times_(times), commTime_(mpi::clone(commTime)),
    subWinLength_(),
    nmembers_(members.size()), localmembers_(), members_(members), commEns_(mpi::clone(commEns)),
    mymembers_(),
    dataset_()
{
  localtimes_ = ntimes_ / commTime_.size();
  localmembers_ = nmembers_ / commEns_.size();
  for (size_t jj = 0; jj < localmembers_; ++jj) {
    size_t indx = localmembers_ * commEns_.rank() + jj;
    ASSERT(indx < members_.size());
    mymembers_.push_back(members_[indx]);
  }
  if (ntimes_ > 1) {
    subWinLength_ = times_[1] - times_[0];
    for (size_t jt = 2; jt < ntimes_; ++jt) {
      ASSERT(times_[jt] - times_[jt-1] == subWinLength_);
    }
  }
}

// -----------------------------------------------------------------------------

template<typename DATA>
DataSetBase<DATA>::DataSetBase(const eckit::mpi::Comm & commTime,
                               const eckit::mpi::Comm & commEns)
  : ntimes_(), localtimes_(), times_(), commTime_(mpi::clone(commTime)),
    subWinLength_(),
    nmembers_(), localmembers_(), members_(), commEns_(mpi::clone(commEns)),
    mymembers_(),
    dataset_()
{}

// -----------------------------------------------------------------------------

template <typename DATA>
DataSetBase<DATA> & DataSetBase<DATA>::operator=(const DataSetBase<DATA> & other) {
  this->check_consistency(other);
  for (size_t jj = 0; jj < dataset_.size(); ++jj) dataset_[jj] = other.dataset_[jj];
  return *this;
}

// -----------------------------------------------------------------------------

template <typename DATA>
void DataSetBase<DATA>::check_consistency() const {
  bool error = false;
  if (ntimes_ <= 0) {
    Log::error() << classname() << " error ntimes = " << ntimes_ << std::endl;
    error = true;
  }
  if (nmembers_ <= 0) {
    Log::error() << classname() << " error nmembers = " << nmembers_ << std::endl;
    error = true;
  }
  if (localtimes_ <= 0) {
    Log::error() << classname() << " error localtimes = " << localtimes_ << std::endl;
    error = true;
  }
  if (localmembers_ <= 0) {
    Log::error() << classname() << " error localmembers_ = " << localmembers_ << std::endl;
    error = true;
  }
  if (dataset_.size() != localtimes_ * localmembers_) {
    Log::error() << classname() << " error size : " << dataset_.size() << " " << localtimes_
                 << " " << localmembers_ << std::endl;
    error = true;
  }
  if (nmembers_ % commEns_.size() != 0) {
    Log::error() << classname() << " error ens comm : " << nmembers_ << " " << commEns_.size()
                 << std::endl;
    error = true;
  }
  if (ntimes_ % commTime_.size() != 0) {
    Log::error() << classname() << " error comm time : " << ntimes_ << " " << commTime_.size()
                 << std::endl;
    error = true;
  }
  size_t indx = localtimes_;
  for (size_t jm = 1; jm < localmembers_; ++jm) {
    for (size_t jt = 0; jt < localtimes_; ++jt) {
      if (dataset_.at(jt).validTime() != dataset_.at(indx).validTime()) {
        Log::error() << classname() << " error times : " << jt << " " << indx
                     << " " << dataset_.at(jt).validTime()
                     << " " << dataset_.at(indx).validTime() << std::endl;
        error = true;
      }
      ++indx;
    }
  }
  for (size_t jj = 1; jj < dataset_.size(); ++jj) {
    if (dataset_[0].variables() != dataset_[jj].variables()) {
      Log::error() << classname() << " error variables : " << jj << " "
                   << dataset_[0].variables() << dataset_[jj].variables() << std::endl;
      error = true;
    }
  }
  if (error) throw eckit::BadValue("Inconsistent data set", Here());
}

// -----------------------------------------------------------------------------

template <typename DATA>
void DataSetBase<DATA>::check_consistency(const DataSetBase<DATA> & other) const {
  this->check_consistency();
  other.check_consistency();

  bool error = false;
  if (ntimes_ != other.ntimes_) {
    Log::error() << classname() << " inconsistent ntimes : " << ntimes_
                 << " " << other.ntimes_ << std::endl;
    error = true;
  }
  if (nmembers_ != other.nmembers_) {
    Log::error() << classname() << " inconsistent nmembers : " << nmembers_
                 << " " << other.nmembers_ << std::endl;
    error = true;
  }
  if (localtimes_ != other.localtimes_) {
    Log::error() << classname() << " inconsistent localtimes : " << localtimes_
                 << " " << other.localtimes_ << std::endl;
    error = true;
  }
  if (localmembers_ != other.localmembers_) {
    Log::error() << classname() << " inconsistent localmembers : " << localmembers_
                 << " " << other.localmembers_ << std::endl;
    error = true;
  }
  if (mymembers_ != other.mymembers_) {
    Log::error() << classname() << " inconsistent mymembers : " << mymembers_
                 << " " << other.mymembers_ << std::endl;
    error = true;
  }
  for (size_t jt = 0; jt < localtimes_; ++jt) {
    if (dataset_.at(jt).validTime() != other.dataset_.at(jt).validTime()) {
    Log::error() << classname() << " inconsistent times : " << jt
                 << " " << dataset_.at(jt).validTime()
                 << " " << other.dataset_.at(jt).validTime() << std::endl;
    error = true;
  }
  }
  if (commTime_.size() != other.commTime_.size()) {
    Log::error() << classname() << " inconsistent comm time size: " << commTime_.size()
                 << " " << other.commTime_.size() << std::endl;
    error = true;
  }
  if (commTime_.rank() != other.commTime_.rank()) {
    Log::error() << classname() << " inconsistent comm time rank: " << commTime_.rank()
                 << " " << other.commTime_.rank() << std::endl;
    error = true;
  }
  if (commEns_.size() != other.commEns_.size()) {
    Log::error() << classname() << " inconsistent comm ens size: " << commEns_.size()
                 << " " << other.commEns_.size() << std::endl;
    error = true;
  }
  if (commEns_.rank() != other.commEns_.rank()) {
    Log::error() << classname() << " inconsistent comm ens rank: " << commEns_.rank()
                 << " " << other.commEns_.rank() << std::endl;
    error = true;
  }
  if (error) {
    throw eckit::BadValue("Inconsistent data sets", Here());
  }
}

// -----------------------------------------------------------------------------

template<typename DATA>
const std::vector<util::DateTime> DataSetBase<DATA>::validTimes() const {
  std::vector<util::DateTime> times(localtimes_);
  for (size_t jt = 0; jt < localtimes_; ++jt) {
    times[jt] = dataset_.at(jt).validTime();
  }
  return times;
}

// -----------------------------------------------------------------------------

template<typename DATA>
void DataSetBase<DATA>::shift_forward() {
  Log::trace() << "DataSetBase::shift_forward start" << std::endl;
  static int tag = 15357;
  size_t mytime = commTime_.rank();

// Send values at end of my subwindow to next subwindow
  if (mytime + 1 < commTime_.size()) {
    util::Timer timer("oops::mpi", "send");
    std::vector<double> sendbuf;
    for (size_t jm = 0; jm < localmembers_; ++jm) {
      data(localtimes_ - 1, jm).serialize(sendbuf);
    }
    commTime_.send(sendbuf.data(), sendbuf.size(), mytime+1, tag);
  }

// Shift local data
  for (size_t jm = 0; jm < localmembers_; ++jm) {
    for (size_t jt = 1; jt < localtimes_; ++jt) {
      data(jt, jm) = data(jt - 1, jm);
    }
  }

// Receive values at beginning of my subwindow from previous subwindow
  if (mytime > 0) {
    util::Timer timer("oops::mpi", "receive");
    size_t sz = localmembers_ * data(localtimes_ - 1, 0).serialSize();
    std::vector<double> recvbuf(sz);
    eckit::mpi::Status status = commTime_.receive(recvbuf.data(), sz, mytime-1, tag);
    size_t ii = 0;
    for (size_t jm = 0; jm < localmembers_; ++jm) {
      data(0, jm).deserialize(recvbuf, ii);
    }
    ASSERT(ii == sz);
  } else {
    for (size_t jm = 0; jm < localmembers_; ++jm) {
      data(0, jm).zero(times_[0]);
    }
  }

  ++tag;
  Log::trace() << "DataSetBase::shift_forward done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename DATA>
void DataSetBase<DATA>::shift_backward() {
  Log::trace() << "DataSetBase::shift_backward start" << std::endl;
  static int tag = 30951;
  size_t mytime = commTime_.rank();

// Send values at start of my subwindow to previous subwindow
  if (mytime > 0) {
    util::Timer timer("oops::mpi", "send");
    std::vector<double> sendbuf;
    for (size_t jm = 0; jm < localmembers_; ++jm) {
      data(0, jm).serialize(sendbuf);
    }
    commTime_.send(sendbuf.data(), sendbuf.size(), mytime-1, tag);
  }

// Shift local data
  for (size_t jm = 0; jm < localmembers_; ++jm) {
    for (size_t jt = 0; jt < localtimes_ - 1; ++jt) {
      data(jt, jm) = data(jt + 1, jm);
    }
  }

// Receive values at end of my subwindow from next subwindow
  if (mytime + 1 < commTime_.size()) {
    util::Timer timer("oops::mpi", "receive");
    size_t sz = localmembers_ * data(localtimes_ - 1, 0).serialSize();
    std::vector<double> recvbuf(sz);
    eckit::mpi::Status status = commTime_.receive(recvbuf.data(), sz, mytime+1, tag);
    size_t ii = 0;
    for (size_t jm = 0; jm < localmembers_; ++jm) {
      data(localtimes_ - 1, jm).deserialize(recvbuf, ii);
    }
    ASSERT(ii == sz);
  } else {
    for (size_t jm = 0; jm < localmembers_; ++jm) {
      data(localtimes_ - 1, jm).zero();
      data(localtimes_ - 1, jm).updateTime(subWinLength_);
    }
  }

  ++tag;
  Log::trace() << "DataSetBase::shift_backward done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename DATA>
void DataSetBase<DATA>::print(std::ostream & os) const {
//  os << classname() << ": " << ntimes_ << " times, " << localtimes_ << " local." << std::endl;
//  os << classname() << ": " << nmembers_ << " members, " << localmembers_ << " local."
//                    << std::endl;
//  os << classname() << ": times: " << dataset_.at(0).validTime();
//  for (size_t jt = 1; jt < localtimes_; ++jt) os << dataset_.at(jt).validTime();
//  os << std::endl;
//  os << classname() << ": members: " << mymembers_ << std::endl;
  if (commTime_.size() > 1) {
    if (dataset_.size() > 1) throw eckit::NotImplemented("gatherprint not good enough", Here());
    gatherPrint(os, dataset_[0], commTime_);
  } else {
    for (const DATA & data : dataset_) os << data;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
