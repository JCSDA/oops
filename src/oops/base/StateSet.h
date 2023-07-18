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

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/DataSetBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template<typename MODEL> class StateSet : public DataSetBase<State<MODEL>> {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef DataSetBase<State_>        Base_;

 public:
  StateSet(const std::vector<util::DateTime> &, const eckit::mpi::Comm &,
           const std::vector<int> & ens = {0},
           const eckit::mpi::Comm & commEns = oops::mpi::myself());
  StateSet(const Geometry_ &, const eckit::Configuration &,
           const eckit::mpi::Comm & commTime = oops::mpi::myself(),
           const eckit::mpi::Comm & commEns = oops::mpi::myself());
  StateSet(const Geometry_ &, const StateSet &);
  StateSet(const StateSet &) = default;
  virtual ~StateSet() = default;

  const Geometry_ & geometry() const {return Base_::dataset_[0].geometry();}

  void read(const Geometry_ &, const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

 protected:
  const std::vector<State_> & states() const {return Base_::dataset_;}
  std::vector<State_> & states() {return Base_::dataset_;}

 private:
  std::string classname() const {return "StateSet";}
};

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const std::vector<util::DateTime> & times,
                          const eckit::mpi::Comm & commTime,
                          const std::vector<int> & ens,
                          const eckit::mpi::Comm & commEns)
  : DataSetBase<State_>(times, commTime, ens, commEns)
{
  Log::trace() << "StateSet::StateSet" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol, const eckit::Configuration & config,
                          const eckit::mpi::Comm & commTime, const eckit::mpi::Comm & commEns)
  : DataSetBase<State_>(commTime, commEns)
{
  Log::trace() << "StateSet::StateSet start " << config << std::endl;
  this->read(resol, config);
  Log::trace() << "StateSet::StateSet done" << *this << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol, const StateSet & other)
  : DataSetBase<State_>(other.times_, other.commTime_, other.members_, other.commEns_)
{
  Log::trace() << "StateSet::StateSet chres start" << std::endl;
  for (size_t jj = 0; jj < other.size(); ++jj) {
    Base_::dataset_.emplace_back(State_(resol, other[jj]));
  }
  Log::trace() << "StateSet::StateSet chres done" << *this << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StateSet<MODEL>::read(const Geometry_ & resol, const eckit::Configuration & config) {
  Log::trace() << "StateSet::read start " << config << std::endl;
  ASSERT(Base_::dataset_.size() == 0);

  if (config.has("members from template")) {
    throw eckit::BadParameter("Members template should be expanded", Here());
  }

  std::vector<eckit::LocalConfiguration> ensconfs;
  if (config.has("members")) {
    ensconfs = config.getSubConfigurations("members");
  } else {
    ASSERT(Base_::commEns_.size() == 1);
    ensconfs.emplace_back(eckit::LocalConfiguration(config));
  }
  Base_::nmembers_ = ensconfs.size();
  Base_::localmembers_ = Base_::nmembers_ / Base_::commEns_.size();
  Base_::mymembers_.clear();
  for (size_t jm = 0; jm < Base_::localmembers_; ++jm) {
    Base_::mymembers_.push_back(Base_::commEns_.rank() * Base_::localmembers_ + jm);
  }

  for (size_t jm = 0; jm < Base_::localmembers_; ++jm) {
    const eckit::LocalConfiguration & mconf = ensconfs.at(Base_::mymembers_[jm]);
    std::vector<eckit::LocalConfiguration> confs;
    if (mconf.has("states")) {
      confs = mconf.getSubConfigurations("states");
    } else {
      confs = {mconf};
    }
    Base_::times_.clear();
    Base_::ntimes_ = confs.size();
    Base_::localtimes_ = Base_::ntimes_ / Base_::commTime_.size();
    for (size_t jt = 0; jt < Base_::localtimes_; ++jt) {
      const size_t it = Base_::commTime_.rank() * Base_::localtimes_ + jt;
      Base_::dataset_.emplace_back(State_(resol, confs.at(it)));
      Base_::times_.push_back(Base_::dataset_.back().validTime());
      ASSERT(jm * Base_::localtimes_ + jt == Base_::dataset_.size() - 1);
    }
    ASSERT(Base_::times_.size() == Base_::localtimes_);
  }

  mpi::allGatherv(Base_::commTime_, Base_::times_);
  ASSERT(Base_::times_.size() == Base_::ntimes_);

  this->check_consistency();

  Log::trace() << "StateSet::read done" << *this << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StateSet<MODEL>::write(const eckit::Configuration & config) const {
  if (Base_::nmembers_ > 1) throw eckit::NotImplemented("Ensemble write not implented", Here());

  if (config.has("states")) {
    std::vector<eckit::LocalConfiguration> confs;
    config.get("states", confs);
    ASSERT(confs.size() == Base_::ntimes_);
    for (size_t jt = 0; jt < Base_::localtimes_; ++jt) {
      size_t it = Base_::commTime_.rank() * Base_::localtimes_ + jt;
      if (config.has("member")) confs[it].set("member", config.getInt("member"));
      Base_::dataset_[jt].write(confs[it]);
    }
  } else {
    ASSERT(Base_::dataset_.size() == 1);
    Base_::dataset_[0].write(config);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
