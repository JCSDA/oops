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
  typedef Geometry<MODEL>     Geometry_;
  typedef Increment<MODEL>    Increment_;
  typedef StateSet<MODEL>     States_;

 public:
  IncrementSet(const Geometry_ &, const Variables &, const std::vector<util::DateTime> &,
               const eckit::mpi::Comm &,
               const std::vector<int> & ens = {0},
               const eckit::mpi::Comm & commEns = oops::mpi::myself());
  IncrementSet(const IncrementSet &, const bool copy = true);
  IncrementSet(const Geometry_ &, const IncrementSet &);
  IncrementSet(const Geometry_ &, const Variables &, const std::vector<util::DateTime> &,
               const eckit::Configuration &,
               const eckit::mpi::Comm & commTime = oops::mpi::myself(),
               const eckit::mpi::Comm & commEns = oops::mpi::myself());
  IncrementSet(const Geometry_ &, const Variables &, const States_ &);
  IncrementSet(const Geometry_ &, const Variables &, States_ &, const bool clearStates = false);
  virtual ~IncrementSet() = default;

  void zero();
  void random();
  IncrementSet & operator+=(const IncrementSet &);
  IncrementSet & operator-=(const IncrementSet &);
  IncrementSet & operator*=(const double);
  void axpy(const double, const IncrementSet &);
  double dot_product_with(const IncrementSet &) const;
  void schur_product_with(const IncrementSet &);

  void diff(const States_ &, const States_ &);

  IncrementSet ens_mean() const;
  IncrementSet ens_stddev() const;

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
      this->dataset().emplace_back(new Increment_(resol, vars, times[mytime + jt]));
    }
  }

  this->check_consistency();

  Log::trace() << "IncrementSet::IncrementSet done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementSet<MODEL>::IncrementSet(const IncrementSet & other, const bool copy)
  : IncrementSet(other.geometry(), other.variables(), other.times(),
                 other.commTime(), other.members(), other.commEns())
{
  Log::trace() << "IncrementSet::IncrementSet copy" << std::endl;
  if (copy) {
    for (size_t jj = 0; jj < other.size(); ++jj) {
      (*this)[jj] = other[jj];
    }
  }
  this->check_consistency();
  Log::trace() << "IncrementSet::IncrementSet copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementSet<MODEL>::IncrementSet(const Geometry_ & resol, const IncrementSet & other)
  : IncrementSet(resol, other.variables(), other.times(), other.commTime(),
                 other.members(), other.commEns())
{
  Log::trace() << "IncrementSet::IncrementSet chres" << std::endl;
  for (size_t jj = 0; jj < other.size(); ++jj) {
    (*this)[jj] = Increment_(resol, other[jj]);
  }
  this->check_consistency();
  Log::trace() << "IncrementSet::IncrementSet chres done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementSet<MODEL>::IncrementSet(const Geometry_ & resol, const Variables & vars,
                                  const std::vector<util::DateTime> & times,
                                  const eckit::Configuration & config,
                                  const eckit::mpi::Comm & commTime,
                                  const eckit::mpi::Comm & commEns)
  : DataSetBase<Increment_, Geometry_>(commTime, commEns)
{
  Log::trace() << "IncrementSet::IncrementSet read start " << config << std::endl;

  std::vector<eckit::LocalConfiguration> locals = this->configure(config);

  size_t mytime = this->local_time_size() * commTime.rank();
  size_t indx = 0;
  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      this->dataset().emplace_back(std::make_unique<Increment_>(resol,
                                                    vars, times[mytime + jt]));
      this->dataset().back()->read(locals.at(indx));
      ++indx;
    }
  }

  this->sync_times();
  this->check_consistency();

  Log::trace() << "IncrementSet::IncrementSet read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementSet<MODEL>::IncrementSet(const Geometry_ & resol, const Variables & vars,
                                  const States_ & xx)
  : DataSetBase<Increment_, Geometry_>(xx.times(), xx.commTime(), xx.members(), xx.commEns())
{
  Log::trace() << "IncrementSet::IncrementSet from States" << std::endl;
  size_t mytime = this->local_time_size() * this->commTime().rank();
  const std::vector<util::DateTime> times = xx.times();
  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      const util::DateTime time = times[mytime + jt];
      this->dataset().emplace_back(new Increment_(resol, vars, time));
      (*this)(jt, jm).transfer_from_state(xx(jt, jm));
    }
  }

  this->check_consistency();

  Log::trace() << "IncrementSet::IncrementSet from States done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementSet<MODEL>::IncrementSet(const Geometry_ & resol, const Variables & vars,
                                  States_ & xx, const bool clearStates)
  : DataSetBase<Increment_, Geometry_>(xx.times(), xx.commTime(), xx.members(), xx.commEns())
{
  Log::trace() << "IncrementSet::IncrementSet from States" << std::endl;
  size_t mytime = this->local_time_size() * this->commTime().rank();
  const std::vector<util::DateTime> times = xx.times();
  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      const util::DateTime time = times[mytime + jt];
      this->dataset().emplace_back(new Increment_(resol, vars, time));
      (*this)(jt, jm).transfer_from_state(xx(jt, jm));
      if (clearStates) xx.clear(jt, jm);
    }
  }
  if (clearStates) xx.clear();

  this->check_consistency();

  Log::trace() << "IncrementSet::IncrementSet from States done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementSet<MODEL>::diff(const States_ & xx1, const States_ & xx2) {
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
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj].zero();
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void IncrementSet<MODEL>::random() {
  this->check_consistency();
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj].random();
  }
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
  this->check_consistency(other, false);
  ASSERT(other.ens_size() == this->ens_size() || other.ens_size() == 1);

  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    size_t im = jm;
    if (other.ens_size() == 1) im = 0;
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      (*this)(jt, jm) += other(jt, im);
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
IncrementSet<MODEL> & IncrementSet<MODEL>::operator-=(const IncrementSet<MODEL> & other) {
  this->check_consistency(other, false);
  ASSERT(other.ens_size() == this->ens_size() || other.ens_size() == 1);

  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    size_t im = jm;
    if (other.ens_size() == 1) im = 0;
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      (*this)(jt, jm) -= other(jt, im);
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
IncrementSet<MODEL> & IncrementSet<MODEL>::operator*=(const double zz) {
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj] *= zz;
  }
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

template<typename MODEL>
IncrementSet<MODEL> IncrementSet<MODEL>::ens_mean() const {
  Log::trace() << "IncrementSet::ens_mean start" << std::endl;
  if (this->commEns().size() > 1) {
    throw eckit::NotImplemented("IncrementSet::ens_mean not implemented for distributed ensembles",
                                Here());
  }

  IncrementSet<MODEL> mean(this->geometry(), this->variables(), this->times(), this->commTime());
  const double fact = 1.0 / static_cast<double>(this->ens_size());
  for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
    for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
      mean[jt] += (*this)(jt, jm);
    }
    mean[jt] *= fact;
  }

  Log::trace() << "IncrementSet::ens_mean done" << std::endl;
  return mean;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementSet<MODEL> IncrementSet<MODEL>::ens_stddev() const {
  Log::trace() << "IncrementSet::ens_stddev start" << std::endl;
  ASSERT(this->ens_size() > 1);
  if (this->commEns().size() > 1) {
    throw eckit::NotImplemented(
      "IncrementSet::ens_stddev not implemented for distributed ensembles",
      Here());
  }

  IncrementSet stddev(this->geometry(), this->variables(), this->times(), this->commTime());
  // Ensemble mean
  IncrementSet ensmean = this->ens_mean();

  // Calculate the variance
  const double rr = 1.0/(static_cast<double>(this->ens_size()) - 1.0);
  for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
    for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    Increment_ pert((*this)(jt, jm));
    pert -= ensmean[jt];
    pert.schur_product_with(pert);
    stddev[jt].axpy(rr, pert);
    }
    // Square root to get standard deviation
    stddev[jt].fieldSet().sqrt();
    stddev[jt].synchronizeFields();
  }

  Log::trace() << "StateSet::stddev done" << std::endl;
  return stddev;
}

// -----------------------------------------------------------------------------
}  // namespace oops
