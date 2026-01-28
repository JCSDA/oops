/*
 * (C) Copyright 2023 UCAR
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/DataSetBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template<typename MODEL>
class StateSet : public DataSetBase< State<MODEL>, Geometry<MODEL> > {
  typedef Geometry<MODEL>                 Geometry_;
  typedef Increment<MODEL>                Increment_;
  typedef State<MODEL>                    State_;

 public:
  StateSet(const Geometry_ &, const Variables &,
           const std::vector<util::DateTime> &, const eckit::mpi::Comm &,
           const std::vector<int> & ens = {0},
           const eckit::mpi::Comm & commEns = oops::mpi::myself());
  StateSet(const Geometry_ &, const eckit::Configuration &,
           const eckit::mpi::Comm & commTime = oops::mpi::myself(),
           const eckit::mpi::Comm & commEns = oops::mpi::myself());
  // create a StateSet variable from a std::vector of State variables distributed
  // across communicators
  StateSet(const Geometry_ &, const StateSet &);
  StateSet(const StateSet &) = default;
  StateSet(const std::vector<StateSet> &, const int &);
  StateSet(const std::vector<util::DateTime> & times,
           const eckit::mpi::Comm & commTime,
           const std::vector<int> & members, const eckit::mpi::Comm & commEns,
           const std::vector<std::shared_ptr<State_>> & dataset);

  // Calculate the ensemble mean and return a new StateSet variable
  StateSet ens_mean() const;
  // Collect distributed states and return a local subset
  std::vector<StateSet> transpose(const eckit::mpi::Comm & global,
           const Geometry_ & DAgeometry, const int ensNum) const;
  /// Zero
  void zero();
  /// Accumulator
  void accumul(const double &, const StateSet &);
  virtual ~StateSet() = default;

 private:
  std::string classname() const {return "StateSet";}
};

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol,
                          const Variables & vars,
                          const std::vector<util::DateTime> & times,
                          const eckit::mpi::Comm & commTime,
                          const std::vector<int> & ens,
                          const eckit::mpi::Comm & commEns)
  : DataSetBase<State_, Geometry_>(times, commTime, ens, commEns)
{
  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      this->dataset().emplace_back(new State_(resol, vars, times[jt]));
    }
  }
  this->sync_times();
  this->check_consistency();
  Log::trace() << "StateSet::StateSet" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol, const eckit::Configuration & config,
                          const eckit::mpi::Comm & commTime, const eckit::mpi::Comm & commEns)
  : DataSetBase<State_, Geometry_>(commTime, commEns)
{
  Log::trace() << "StateSet::StateSet read start " << config << std::endl;

// get vector of local configurations
  std::vector<eckit::LocalConfiguration> locals = this->configure(config);

  size_t indx = 0;
  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      this->dataset().emplace_back(new State_(resol, locals.at(indx)));
      ++indx;
    }
  }

  this->sync_times();
  this->check_consistency();

  Log::trace() << "StateSet::StateSet read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol, const StateSet & other)
  : DataSetBase<State_, Geometry_>(other.times(), other.commTime(),
                                   other.members(), other.commEns())
{
  Log::trace() << "StateSet::StateSet chres start" << std::endl;
  for (size_t jj = 0; jj < other.size(); ++jj) {
    this->dataset().emplace_back(std::make_unique<State_>(resol, other[jj]));
  }
  Log::trace() << "StateSet::StateSet chres done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
std::vector<StateSet<MODEL> > StateSet<MODEL>::transpose(const eckit::mpi::Comm & global,
           const Geometry_ & DAgeometry, const int ensNum) const
{
/* This method collects parts of the distributed StateSet and places all ensemble
   member states in a smaller patch (1/N the size of Forecast geometry) of a StateSet 
   held in the local_ensemble. It is essentially a transpose of a distributed StateSet
   to a locally held vector of StateSets. The std::vector of StateSets is used here because
   the LocalEnsemble infrastructure still expects that rather than a normal, single StateSet 
   variable. If that infrastructure changes, the localize call below will support the new 
   approach and this should be deprecated. The DAgeometry should be have a decomposition 
   that is spread across N (number of ensemble members) times the number of MPI tasks that 
   the forecast geometry decomposition. In other words, if the forecast geometry has a 
   layout of [4,4] and there are 9 ensemble members, the DA geometry should have a 
   layout that multiplies to 4*4*9 or something like 12,12. Note that the resolution
   of both geometries is the same (e.g. C48, C96, etc.). Just the decomposition
   is different between the geometries.
*/
  std::vector<StateSet<MODEL> > local;
  std::vector<int> local_ens;
  local_ens.push_back(1);


  /* transpose stateSet to get all ensemble members on a 1/N size patch of geometry */
  for (size_t jm = 0; jm < this->ens_size(); ++jm) {
    local.emplace_back(StateSet(DAgeometry, this->variables(),
       this->times(), this->commTime(), local_ens, oops::mpi::myself()));
    local[jm](0, 0).transpose((*this)(0, 0).state(), global, ensNum, jm);
    local[jm].sync_times();
  }
  return(local);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const std::vector<util::DateTime> & times,
                          const eckit::mpi::Comm & commTime,
                          const std::vector<int> & members, const eckit::mpi::Comm & commEns,
                          const std::vector<std::shared_ptr<State_>> & dataset)
  : DataSetBase<State_, Geometry_>(times, commTime, members, commEns, dataset) {
  Log::trace() << "StateSet::StateSet from dataset done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL> StateSet<MODEL>::ens_mean() const {
  Log::trace() << "StateSet::ens_mean start" << std::endl;
  StateSet<MODEL> mean(this->geometry(), this->variables(), this->times(), this->commTime());
  const double fact = 1.0 / static_cast<double>(this->ens_size());
  for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
    for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
      // get sum of values on local ensemble
      mean[jt].accumul(fact, (*this)(jt, jm));
    }
    if (this->commEns().size() > 1) {
      // sum up values across ensemble communicator
      oops::mpi::allReduceInPlace(this->commEns(), mean[jt].fieldSet().fieldSet());
      mean[jt].synchronizeFields();
    }
  }
  Log::trace() << "StateSet::ens_mean done" << std::endl;
  return mean;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StateSet<MODEL>::zero() {
  Log::trace() << "StateSet<MODEL>::zero starting" << std::endl;
  for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
    for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
      (*this)(jt, jm).zero();
    }
  }
  Log::trace() << "StateSet<MODEL>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------

template<typename MODEL>
void StateSet<MODEL>::accumul(const double & zz, const StateSet & xx) {
  Log::trace() << "StateSet<MODEL>::accumul starting" << std::endl;
  for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
    for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
      (*this)(jt, jm).accumul(zz, xx(jt, jm));
    }
  }
  Log::trace() << "StateSet<MODEL>::accumul done" << std::endl;
}

}  // namespace oops
