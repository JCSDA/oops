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

template<typename MODEL>
class StateSet : public DataSetBase< State<MODEL>, Geometry<MODEL> > {
  typedef Geometry<MODEL>                 Geometry_;
  typedef State<MODEL>                    State_;

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

 protected:
  const std::vector<State_> & states() const {return this->dataset();}
  std::vector<State_> & states() {return this->dataset();}

 private:
  std::string classname() const {return "StateSet";}
};

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const std::vector<util::DateTime> & times,
                          const eckit::mpi::Comm & commTime,
                          const std::vector<int> & ens,
                          const eckit::mpi::Comm & commEns)
  : DataSetBase<State_, Geometry_>(times, commTime, ens, commEns)
{
  Log::trace() << "StateSet::StateSet" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol, const eckit::Configuration & config,
                          const eckit::mpi::Comm & commTime, const eckit::mpi::Comm & commEns)
  : DataSetBase<State_, Geometry_>(commTime, commEns)
{
  Log::trace() << "StateSet::StateSet read start " << config << std::endl;
  this->read(resol, config);
  Log::trace() << "StateSet::StateSet read done" << *this << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol, const StateSet & other)
  : DataSetBase<State_, Geometry_>(other.times(), other.commTime(),
                                   other.members(), other.commEns())
{
  Log::trace() << "StateSet::StateSet chres start" << std::endl;
  for (size_t jj = 0; jj < other.size(); ++jj) {
    this->dataset().emplace_back(State_(resol, other[jj]));
  }
  Log::trace() << "StateSet::StateSet chres done" << *this << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
