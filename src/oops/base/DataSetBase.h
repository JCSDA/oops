/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/mpi/mpi.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/gatherPrint.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

class Variables;

// -----------------------------------------------------------------------------

template<typename DATA, typename GEOM>
class DataSetBase : public util::Printable {
 public:
  virtual ~DataSetBase() = default;

  DataSetBase & operator=(const DataSetBase &);

  bool is_3d() const {return nmembers_ == 1 && ntimes_ == 1;}
  bool is_4d() const {return nmembers_ == 1;}

  size_t ens_size() const {return nmembers_;}
  const std::vector<int> & members() const {return allmembers_;}

  size_t time_size() const {return ntimes_;}
  const std::vector<util::DateTime> & times() const {
    sync_times();
    return alltimes_;
  }

  const size_t & local_ens_size() const {return localmembers_;}
  const size_t & local_time_size() const {return localtimes_;}
  DATA & operator()(const size_t it, const size_t im) {return this->data(it, im);}
  const DATA & operator()(const size_t it, const size_t im) const {return this->data(it, im);}

  size_t size() const {return dataset_.size();}
  DATA & operator[](const int ii) {return *dataset_[ii];}
  const DATA & operator[](const int ii) const {return *dataset_[ii];}

  const GEOM & geometry() const {return dataset_[0]->geometry();}
  const Variables & variables() const {this->check_consistency(); return dataset_[0]->variables();}
  const std::vector<util::DateTime> validTimes() const;

  const eckit::mpi::Comm & commEns() const {return commEns_;}
  const eckit::mpi::Comm & commTime() const {return commTime_;}

  void read(const GEOM &, const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

  void shift_forward();
  void shift_backward();

  void clear(const size_t it, const size_t im) {dataset_.at(im * localtimes_ + it).reset();}
  void clear();

 protected:
  DataSetBase(const std::vector<util::DateTime> &, const eckit::mpi::Comm &,
              const std::vector<int> &, const eckit::mpi::Comm &);
  DataSetBase(const eckit::mpi::Comm &, const eckit::mpi::Comm &);
  DataSetBase(const DataSetBase &);
  std::vector<eckit::LocalConfiguration> configure(const eckit::Configuration &);
  void sync_times() const;

  void check_consistency() const;
  void check_consistency(const DataSetBase &, const bool strict_members = true) const;

  std::vector<std::unique_ptr<DATA>> & dataset() {return dataset_;}
  const std::vector<std::unique_ptr<DATA>> & dataset() const {return dataset_;}

 private:
  DATA & data(const size_t it, const size_t im) {return *dataset_.at(im * localtimes_ + it);}
  const DATA & data(const size_t it, const size_t im) const
                                                {return *dataset_.at(im * localtimes_ + it);}

  size_t ntimes_;
  size_t localtimes_;
  mutable std::vector<util::DateTime> alltimes_;
  const eckit::mpi::Comm & commTime_;
  util::Duration subWinLength_;

  size_t nmembers_;
  size_t localmembers_;
  std::vector<int> allmembers_;
  const eckit::mpi::Comm & commEns_;
  std::vector<size_t> mymembers_;

  std::vector<std::unique_ptr<DATA>> dataset_;

  void print(std::ostream &) const;
  virtual std::string classname() const = 0;
};

// -----------------------------------------------------------------------------

template<typename DATA, typename GEOM>
DataSetBase<DATA, GEOM>::DataSetBase(const std::vector<util::DateTime> & times,
                                     const eckit::mpi::Comm & commTime,
                                     const std::vector<int> & members,
                                     const eckit::mpi::Comm & commEns)
  : ntimes_(times.size()), localtimes_(), alltimes_(times), commTime_(mpi::clone(commTime)),
    subWinLength_(),
    nmembers_(members.size()), localmembers_(), allmembers_(members), commEns_(mpi::clone(commEns)),
    mymembers_(),
    dataset_()
{
  localtimes_ = ntimes_ / commTime_.size();
  localmembers_ = nmembers_ / commEns_.size();
  for (size_t jj = 0; jj < localmembers_; ++jj) {
    size_t indx = localmembers_ * commEns_.rank() + jj;
    ASSERT(indx < allmembers_.size());
    mymembers_.push_back(allmembers_[indx]);
  }
  if (ntimes_ > 1) {
    subWinLength_ = alltimes_[1] - alltimes_[0];
    for (size_t jt = 2; jt < ntimes_; ++jt) {
      ASSERT(alltimes_[jt] - alltimes_[jt-1] == subWinLength_);
    }
  }
}

// -----------------------------------------------------------------------------

template<typename DATA, typename GEOM>
DataSetBase<DATA, GEOM>::DataSetBase(const eckit::mpi::Comm & commTime,
                                     const eckit::mpi::Comm & commEns)
  : ntimes_(), localtimes_(), alltimes_(), commTime_(mpi::clone(commTime)),
    subWinLength_(),
    nmembers_(), localmembers_(), allmembers_(), commEns_(mpi::clone(commEns)),
    mymembers_(),
    dataset_()
{}

// -----------------------------------------------------------------------------

template <typename DATA, typename GEOM>
DataSetBase<DATA, GEOM>::DataSetBase(const DataSetBase<DATA, GEOM> & other)
  : ntimes_(other.ntimes_), localtimes_(other.localtimes_), alltimes_(other.alltimes_),
    commTime_(mpi::clone(other.commTime_)), subWinLength_(other.subWinLength_),
    nmembers_(other.nmembers_), localmembers_(other.localmembers_), allmembers_(other.allmembers_),
    commEns_(mpi::clone(other.commEns_)), mymembers_(other.mymembers_),
    dataset_()
{
  for (size_t jj = 0; jj < other.dataset_.size(); ++jj) {
    dataset_.emplace_back(new DATA(*other.dataset_[jj]));
  }
  this->check_consistency(other);
}

// -----------------------------------------------------------------------------

template <typename DATA, typename GEOM>
DataSetBase<DATA, GEOM> & DataSetBase<DATA, GEOM>::operator=(const DataSetBase<DATA, GEOM> & other)
{
  this->check_consistency(other);
  for (size_t jj = 0; jj < dataset_.size(); ++jj) *dataset_[jj] = *other.dataset_[jj];
  return *this;
}

// -----------------------------------------------------------------------------

template <typename DATA, typename GEOM>
void DataSetBase<DATA, GEOM>::check_consistency() const {
  bool error = false;
  if (ntimes_ == 0 && nmembers_ == 0) return;
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
      if (dataset_.at(jt)->validTime() != dataset_.at(indx)->validTime()) {
        Log::error() << classname() << " error times : " << jt << " " << indx
                     << " " << dataset_.at(jt)->validTime()
                     << " " << dataset_.at(indx)->validTime() << std::endl;
        error = true;
      }
      ++indx;
    }
  }
  for (size_t jj = 1; jj < dataset_.size(); ++jj) {
    if (dataset_[0]->variables() != dataset_[jj]->variables()) {
      Log::error() << classname() << " error variables : " << jj << " "
                   << dataset_[0]->variables() << dataset_[jj]->variables() << std::endl;
      error = true;
    }
  }
  if (error) throw eckit::BadValue("Inconsistent data set", Here());
}

// -----------------------------------------------------------------------------

template <typename DATA, typename GEOM>
void DataSetBase<DATA, GEOM>::check_consistency(const DataSetBase<DATA, GEOM> & other,
                                                const bool strict_members) const {
  this->check_consistency();
  other.check_consistency();

  bool error = false;
  if (ntimes_ != other.ntimes_) {
    Log::error() << classname() << " inconsistent ntimes : " << ntimes_
                 << " " << other.ntimes_ << std::endl;
    error = true;
  }
  if (localtimes_ != other.localtimes_) {
    Log::error() << classname() << " inconsistent localtimes : " << localtimes_
                 << " " << other.localtimes_ << std::endl;
    error = true;
  }
  for (size_t jt = 0; jt < localtimes_; ++jt) {
    if (dataset_.at(jt)->validTime() != other.dataset_.at(jt)->validTime()) {
      Log::error() << classname() << " inconsistent times : " << jt
                   << " " << dataset_.at(jt)->validTime()
                   << " " << other.dataset_.at(jt)->validTime() << std::endl;
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
  if (strict_members) {
    if (nmembers_ != other.nmembers_) {
      Log::error() << classname() << " inconsistent nmembers : " << nmembers_
                   << " " << other.nmembers_ << std::endl;
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
  }
  if (error) {
    throw eckit::BadValue("Inconsistent data sets", Here());
  }
}

// -----------------------------------------------------------------------------

template<typename DATA, typename GEOM>
const std::vector<util::DateTime> DataSetBase<DATA, GEOM>::validTimes() const {
  std::vector<util::DateTime> times(localtimes_);
  for (size_t jt = 0; jt < localtimes_; ++jt) {
    times[jt] = dataset_.at(jt)->validTime();
  }
  return times;
}

// -----------------------------------------------------------------------------

template<typename DATA, typename GEOM>
std::vector<eckit::LocalConfiguration>
  DataSetBase<DATA, GEOM>::configure(const eckit::Configuration & config)
{
  Log::trace() << "DataSetBase::configure start " << config << std::endl;

  std::vector<eckit::LocalConfiguration> locals;

  std::vector<eckit::LocalConfiguration> ensconfs;
  if (config.has("members from template")) {
    eckit::LocalConfiguration tmpl(config, "members from template");
    nmembers_ = tmpl.getInt("nmembers");
    const std::string pattern = tmpl.getString("pattern");
    const int zpad = tmpl.getInt("zero padding", 0);
    const std::vector<size_t> except = tmpl.getUnsignedVector("except", {});
    size_t index = tmpl.getUnsigned("start", 1);
    for (size_t jens = 0; jens < nmembers_; ++jens) {
      while (std::count(except.begin(), except.end(), index)) {
        index++;
      }
      eckit::LocalConfiguration conf(tmpl, "template");
      util::seekAndReplace(conf, pattern, index, zpad);
      ensconfs.push_back(conf);
      index++;
    }
  } else if (config.has("members")) {
    ensconfs = config.getSubConfigurations("members");
    nmembers_ = ensconfs.size();
  } else {
    ASSERT(commEns_.size() == 1);
    ensconfs.emplace_back(eckit::LocalConfiguration(config));
    nmembers_ = 1;
  }

  ASSERT(ensconfs.size() == nmembers_);
  allmembers_.resize(nmembers_);
  std::iota(std::begin(allmembers_), std::end(allmembers_), 0);
  localmembers_ = nmembers_ / commEns_.size();
  mymembers_.clear();
  for (size_t jm = 0; jm < localmembers_; ++jm) {
    mymembers_.push_back(commEns_.rank() * localmembers_ + jm);
  }

  for (size_t jm = 0; jm < localmembers_; ++jm) {
    const eckit::LocalConfiguration & mconf = ensconfs.at(mymembers_[jm]);
    std::vector<eckit::LocalConfiguration> confs;
    if (mconf.has("states")) {
      confs = mconf.getSubConfigurations("states");
    } else {
      confs = {mconf};
    }
    if (jm == 0) {
      ntimes_ = confs.size();
      localtimes_ = ntimes_ / commTime_.size();
    } else {
      ASSERT(ntimes_ == confs.size());
      ASSERT(localtimes_ == ntimes_ / commTime_.size());
    }
    for (size_t jt = 0; jt < localtimes_; ++jt) {
      const size_t it = commTime_.rank() * localtimes_ + jt;
      locals.push_back(confs.at(it));
    }
  }

  ASSERT(locals.size() == localtimes_ * localmembers_);

  Log::trace() << "DataSetBase::configure done" << std::endl;
  return locals;
}

// -----------------------------------------------------------------------------

template<typename DATA, typename GEOM>
void DataSetBase<DATA, GEOM>::read(const GEOM & resol, const eckit::Configuration & config) {
  Log::trace() << "DataSetBase::read start " << config << std::endl;

  std::vector<eckit::LocalConfiguration> locals = this->configure(config);

  for (size_t jm = 0; jm < localmembers_; ++jm) {
    for (size_t jt = 0; jt < localtimes_; ++jt) {
      this->data(jt, jm).read(locals.at(localtimes_ * jm + jt));
    }
  }

  this->sync_times();
  this->check_consistency();

  Log::trace() << "DataSetBase::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename DATA, typename GEOM>
void DataSetBase<DATA, GEOM>::write(const eckit::Configuration & config) const {
  Log::trace() << "DataSetBase::write start" << std::endl;
  if (nmembers_ > 1) throw eckit::NotImplemented("Ensemble write not implented", Here());

  if (config.has("states")) {
    std::vector<eckit::LocalConfiguration> confs;
    config.get("states", confs);
    ASSERT(confs.size() == ntimes_);
    for (size_t jt = 0; jt < localtimes_; ++jt) {
      size_t it = commTime_.rank() * localtimes_ + jt;
      if (config.has("member")) confs[it].set("member", config.getInt("member"));
      dataset_[jt]->write(confs[it]);
    }
  } else {
    ASSERT(dataset_.size() == 1);
    dataset_[0]->write(config);
  }
  Log::trace() << "DataSetBase::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename DATA, typename GEOM>
void DataSetBase<DATA, GEOM>::shift_forward() {
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
      data(0, jm).zero(alltimes_[0]);
    }
  }

  ++tag;
  Log::trace() << "DataSetBase::shift_forward done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename DATA, typename GEOM>
void DataSetBase<DATA, GEOM>::shift_backward() {
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

template <typename DATA, typename GEOM>
void DataSetBase<DATA, GEOM>::sync_times() const {
  alltimes_.resize(localtimes_);
  alltimes_ = this->validTimes();
  mpi::allGatherv(commTime_, alltimes_);
  ASSERT(alltimes_.size() == ntimes_);
}

// -----------------------------------------------------------------------------

template <typename DATA, typename GEOM>
void DataSetBase<DATA, GEOM>::clear() {
  alltimes_.clear();
  allmembers_.clear();
  mymembers_.clear();
  dataset_.clear();
  ntimes_ = 0;
  localtimes_ = 0;
  nmembers_ = 0;
  localmembers_ = 0;
}

// -----------------------------------------------------------------------------

template <typename DATA, typename GEOM>
void DataSetBase<DATA, GEOM>::print(std::ostream & os) const {
//  os << classname() << ": " << nmembers_ << " members, " << localmembers_ << " local."
//                    << std::endl;
//  os << classname() << ": members: " << mymembers_ << std::endl;
//  os << classname() << ": " << ntimes_ << " times, " << localtimes_ << " local." << std::endl;
//  os << classname() << ": times: " << alltimes_[0];
//  for (size_t jt = 1; jt < ntimes_; ++jt) os << ", " << alltimes_[jt];
//  os << std::endl;
//  os << classname() << ": valid times: " << dataset_.at(0)->validTime();
//  for (size_t jt = 1; jt < localtimes_; ++jt) os << ", " << dataset_.at(jt)->validTime();
//  os << std::endl;

  if (commTime_.size() > 1 && commEns_.size() > 1)
    throw eckit::NotImplemented("gatherprint not good enough", Here());

  if (commTime_.size() > 1) {
    if (dataset_.size() > 1) throw eckit::NotImplemented("gatherprint not good enough", Here());
    gatherPrint(os, *dataset_[0], commTime_);
  } else if (commEns_.size() > 1) {
    if (dataset_.size() > 1) throw eckit::NotImplemented("gatherprint not good enough", Here());
    gatherPrint(os, *dataset_[0], commEns_);
  } else {
    for (const auto & data : dataset_) os << *data;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
