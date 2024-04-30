/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_INCREMENT4D_H_
#define OOPS_BASE_INCREMENT4D_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------
/// 4D model state Increment (vector of 3D Increments)
template<typename MODEL> class Increment4D : public IncrementSet<MODEL> {
  typedef Geometry<MODEL>    Geometry_;
  typedef Increment<MODEL>   Increment_;

 public:
  /// Constructor for specified times
  Increment4D(const Geometry_ &, const Variables &, const std::vector<util::DateTime> &,
              const eckit::mpi::Comm & commTime = oops::mpi::myself());
  Increment4D(const Increment4D &, const bool copy = true);
  Increment4D(const Geometry_ &, const Increment4D &);

  void ones();
  void dirac(const eckit::Configuration &);

 private:
  std::string classname() const {return "Increment4D";}
};

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment4D<MODEL>::Increment4D(const Geometry_ & geom, const Variables & vars,
                                const std::vector<util::DateTime> & times,
                                const eckit::mpi::Comm & commTime)
  : IncrementSet<MODEL>(geom, vars, times, commTime)
{
  ASSERT(this->is_4d());
  Log::trace() << "Increment4D:Increment4D created." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment4D<MODEL>::Increment4D(const Increment4D & other, const bool copy)
  : IncrementSet<MODEL>(other, copy)
{
  ASSERT(this->is_4d());
  Log::trace() << "Increment4D:Increment4D copy created." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment4D<MODEL>::Increment4D(const Geometry_ & resol, const Increment4D & other)
  : IncrementSet<MODEL>(resol, other)
{
  ASSERT(this->is_4d());
  Log::trace() << "Increment4D:Increment4D created resol." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment4D<MODEL>::ones() {
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj].ones();
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment4D<MODEL>::dirac(const eckit::Configuration & conf) {
  if (this->time_size() == 1) {
    (*this)[0].dirac(conf);
  } else {
    const std::vector<eckit::LocalConfiguration> confs = conf.getSubConfigurations();
    ASSERT(this->time_size() == confs.size());
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      const size_t it = this->commTime().rank() * this->local_time_size() + jt;
      if (confs[it].empty()) {
        (*this)[jt].zero();
      } else {
        (*this)[jt].dirac(confs[it]);
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_INCREMENT4D_H_
