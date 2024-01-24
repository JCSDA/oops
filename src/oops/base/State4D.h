/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// Four dimensional state (vector of 3D States)
template<typename MODEL> class State4D : public StateSet<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;

 public:
  /// The arguments define all states in 4D and their resolution
  State4D(const Geometry_ &, const eckit::Configuration &,
          const eckit::mpi::Comm & commTime = oops::mpi::myself());
  State4D(const Geometry_ &, const State4D &);

  /// Accumulator
  void zero();
  void accumul(const double &, const State4D &);

 private:
  std::string classname() const {return "State4D";}
};


// -----------------------------------------------------------------------------

template<typename MODEL>
State4D<MODEL>::State4D(const Geometry_ & resol, const eckit::Configuration & config,
                        const eckit::mpi::Comm & commTime)
  : StateSet<MODEL>(resol, config, commTime)
{
  ASSERT(this->is_4d());
  Log::trace() << "State4D constructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State4D<MODEL>::State4D(const Geometry_ & resol, const State4D & other)
  : StateSet<MODEL>(resol, other)
{
  ASSERT(this->is_4d());
  Log::trace() << "State4D constructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State4D<MODEL>::zero() {
  Log::trace() << "State4D<MODEL>::zero starting" << std::endl;
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj].zero();
  }
  Log::trace() << "State4D<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State4D<MODEL>::accumul(const double & zz, const State4D & xx) {
  Log::trace() << "State4D<MODEL>::accumul starting" << std::endl;
  ASSERT(xx.size() == this->size());
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj].accumul(zz, xx[jj]);
  }
  Log::trace() << "State4D<MODEL>::accumul done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
