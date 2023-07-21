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

#include "oops/base/DataSetBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/StateSet.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template<typename MODEL>
class IncrementSet : public DataSetBase< Increment<MODEL>, Geometry<MODEL> > {
  typedef Geometry<MODEL>                     Geometry_;
  typedef Increment<MODEL>                    Increment_;

 public:
  IncrementSet(const Geometry_ &, const Variables &, const std::vector<util::DateTime> &,
               const eckit::mpi::Comm &,
               const std::vector<int> & ens = {0},
               const eckit::mpi::Comm & commEns = oops::mpi::myself());
  IncrementSet(const IncrementSet &, const bool copy = true);
  IncrementSet(const Geometry_ &, const IncrementSet &);
  virtual ~IncrementSet() = default;

  void zero();
  void random();
  IncrementSet & operator+=(const IncrementSet &);
  IncrementSet & operator-=(const IncrementSet &);
  IncrementSet & operator*=(const double);
  void axpy(const double, const IncrementSet &);
  double dot_product_with(const IncrementSet &) const;
  void schur_product_with(const IncrementSet &);

  void diff(const StateSet<MODEL> &, const StateSet<MODEL> &);

 protected:
  const std::vector<Increment_> & increments() const {return this->dataset();}
  std::vector<Increment_> & increments() {return this->dataset();}

 private:
  std::string classname() const {return "IncrementSet";}
};

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementSet<MODEL>::IncrementSet(const Geometry_ & resol, const Variables & vars,
                                  const std::vector<util::DateTime> & times,
                                  const eckit::mpi::Comm & commTime,
                                  const std::vector<int> & ens,
                                  const eckit::mpi::Comm & commEns)
  : DataSetBase<Increment_, Geometry_>(times, commTime, ens, commEns)
{
  Log::trace() << "IncrementSet::IncrementSet" << std::endl;

  size_t mytime = this->local_time_size() * commTime.rank();
  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      this->dataset().emplace_back(Increment_(resol, vars, times[mytime + jt]));
    }
  }

  this->check_consistency();

  Log::trace() << "IncrementSet::IncrementSet done" << *this << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementSet<MODEL>::IncrementSet(const IncrementSet & other, const bool copy)
  : IncrementSet(other.geometry(), other.variables(), other.times(),
                 other.commTime(), other.members(), other.commEns())
{
  Log::trace() << "IncrementSet::IncrementSet" << std::endl;
  if (copy) {
    for (size_t jj = 0; jj < other.size(); ++jj) {
      (*this)[jj] = other[jj];
    }
  }
  this->check_consistency();
  Log::trace() << "IncrementSet::IncrementSet copied" << *this << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementSet<MODEL>::IncrementSet(const Geometry_ & resol, const IncrementSet & other)
  : IncrementSet(resol, other.variables(), other.times(), other.commTime(),
                 other.members(), other.commEns())
{
  Log::trace() << "IncrementSet::IncrementSet" << std::endl;
  for (size_t jj = 0; jj < other.size(); ++jj) {
    (*this)[jj] = Increment_(resol, other[jj]);
  }
  this->check_consistency();
  Log::trace() << "IncrementSet::IncrementSet chres done" << *this << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementSet<MODEL>::diff(const StateSet<MODEL> & xx1, const StateSet<MODEL> & xx2) {
//  this->check_consistency(xx1);
//  this->check_consistency(xx2);
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj].diff(xx1[jj], xx2[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void IncrementSet<MODEL>::zero() {
  this->check_consistency();
  for (Increment_ & incr : this->dataset()) incr.zero();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void IncrementSet<MODEL>::random() {
  this->check_consistency();
  for (Increment_ & incr : this->dataset()) incr.random();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
StateSet<MODEL> & operator+=(StateSet<MODEL> & xx, const IncrementSet<MODEL> & dx) {
  Log::trace() << "operator+=(StateSet, IncrementSet) starting" << std::endl;
//  xx.check_consistency(dx);
  for (size_t ii = 0; ii < xx.size(); ++ii) {
    xx[ii] += dx[ii];
  }
  Log::trace() << "operator+=(StateSet, IncrementSet) done" << std::endl;
  return xx;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
IncrementSet<MODEL> & IncrementSet<MODEL>::operator+=(const IncrementSet<MODEL> & other) {
  this->check_consistency(other);
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj] += other[jj];
  }
  return *this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
IncrementSet<MODEL> & IncrementSet<MODEL>::operator-=(const IncrementSet<MODEL> & other) {
  this->check_consistency(other);
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj] -= other[jj];
  }
  return *this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
IncrementSet<MODEL> & IncrementSet<MODEL>::operator*=(const double zz) {
  for (Increment_ & incr : this->dataset()) incr *= zz;
  return *this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void IncrementSet<MODEL>::axpy(const double zz, const IncrementSet<MODEL> & other) {
  this->check_consistency(other);
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj].axpy(zz, other[jj]);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
double IncrementSet<MODEL>::dot_product_with(const IncrementSet<MODEL> & other) const {
  this->check_consistency(other);
  ASSERT(this->is_4d());  // for now
  double zz = 0.0;
  for (size_t jj = 0; jj < this->size(); ++jj) {
    zz += dot_product((*this)[jj], other[jj]);
  }
  this->commTime().allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  this->commEns().allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  return zz;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void IncrementSet<MODEL>::schur_product_with(const IncrementSet<MODEL> & other) {
  this->check_consistency(other);
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj].schur_product_with(other[jj]);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
