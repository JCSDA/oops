/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVATIONS_H_
#define OOPS_BASE_OBSERVATIONS_H_

#include <cstddef>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Departures.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

template<typename MODEL> class Observations;

/// Observations Class.
/*!
 *  Contains observed values or their model equivalents and the associated
 *  observation operator.
 */

// -----------------------------------------------------------------------------
template <typename MODEL> class Observations : public util::Printable {
  typedef Departures<MODEL>          Departures_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  explicit Observations(const ObsSpaces_ &, const std::string name = "");
  Observations(const Observations &);
  Observations(const ObsSpaces_ &, const Observations &);
  ~Observations();
  Observations & operator=(const Observations &);

/// Access
  std::size_t size() const {return obs_.size();}
  ObsVector_ & operator[](const std::size_t ii) {return *obs_.at(ii);}
  const ObsVector_ & operator[](const std::size_t ii) const {return *obs_.at(ii);}

/// Interactions with Departures
  std::vector<std::shared_ptr<ObsVector_> > operator-(const Observations & other) const;
  Observations & operator+=(const Departures_ &);

/// Save observations values
  void save(const std::string &) const;

/// Accumulator
  void zero();
  void accumul(const Observations &);
  Observations & operator*=(const double);

 private:
  void print(std::ostream &) const;

  std::vector<std::shared_ptr<ObsVector_> > obs_;
};

// =============================================================================

template <typename MODEL>
Observations<MODEL>::Observations(const ObsSpaces_ & obsdb,
                                  const std::string name): obs_(0)
{
  for (std::size_t jj = 0; jj < obsdb.size(); ++jj) {
    std::shared_ptr<ObsVector_> tmp(new ObsVector_(obsdb[jj], name));
    obs_.push_back(tmp);
  }
  Log::trace() << "Observations created" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL>::Observations(const Observations & other): obs_(0)
{
// We want deep copy here
  for (std::size_t jj = 0; jj < other.obs_.size(); ++jj) {
    std::shared_ptr<ObsVector_> tmp(new ObsVector_(*other.obs_[jj]));
    obs_.push_back(tmp);
  }
  Log::trace() << "Observations copy-created" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL>::Observations(const ObsSpaces_ & obsdb,
                                  const Observations & other): obs_(0) {
  for (std::size_t jj = 0; jj < other.obs_.size(); ++jj) {
    std::shared_ptr<ObsVector_> tmp(new ObsVector_(obsdb[jj], *other.obs_[jj]));
    obs_.push_back(tmp);
  }
  Log::trace() << "Local observations created" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL>::~Observations() {
  Log::trace() << "Observations destructed" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL> & Observations<MODEL>::operator=(const Observations & rhs) {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    *obs_[jj] = *rhs.obs_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
std::vector<std::shared_ptr<ObsVector<MODEL> > >
Observations<MODEL>::operator-(const Observations & other) const {
  std::vector<std::shared_ptr<ObsVector_> > out;
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    std::shared_ptr<ObsVector_> ovec(new ObsVector_(*obs_[jj]));
    *ovec -= *other.obs_[jj];
    out.push_back(ovec);
  }
  return out;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL> & Observations<MODEL>::operator+=(const Departures_ & dy) {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    *obs_[jj] += dy[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::save(const std::string & name) const {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    obs_[jj]->save(name);
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::zero() {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    obs_[jj]->zero();
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::accumul(const Observations & y) {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    *obs_[jj] += y[jj];
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Observations<MODEL> & Observations<MODEL>::operator *=(const double factor) {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) {
    *obs_[jj] *= factor;
  }
  return *this;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Observations<MODEL>::print(std::ostream & os) const {
  for (std::size_t jj = 0; jj < obs_.size(); ++jj) os << *obs_[jj] << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_OBSERVATIONS_H_
