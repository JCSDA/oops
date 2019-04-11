/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_DEPARTURES_H_
#define OOPS_BASE_DEPARTURES_H_

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/ObsErrorBase.h"
#include "oops/base/ObsOperators.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

template<typename MODEL> class Departures;

/// Difference between two observation vectors.
/*!
 * A departure is the difference between two observations.
 * The archetypal example is \f$ \mathbf{y} - {\cal H}(\mathbf{x}) \f$.
 *
 * Keeping an observation space vector here is necessary for the implementation
 * of generic observation error covariance matrices.
 */

// -----------------------------------------------------------------------------
template <typename MODEL>
class Departures : public util::Printable,
                   public GeneralizedDepartures {
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncr_;
  typedef ObsErrorBase<MODEL>        ObsErrorBase_;
  typedef ObsOperators<MODEL>        ObsOperators_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
// Constructors and destructor
  Departures(const ObsSpaces_ &, const ObsOperators_ &);
  explicit Departures(std::vector<boost::shared_ptr<ObsVector_> >);
  explicit Departures(const Departures &);
  ~Departures();

/// Access
  std::size_t size() const {return dep_.size();}
  ObsVector_ & operator[](const std::size_t ii) {return *dep_.at(ii);}
  const ObsVector_ & operator[](const std::size_t ii) const {return *dep_.at(ii);}

// Linear algebra operators
  Departures & operator=(const Departures &);
  Departures & operator+=(const Departures &);
  Departures & operator-=(const Departures &);
  Departures & operator*=(const double &);
  Departures & operator*=(const Departures &);
  Departures & operator/=(const Departures &);
  void zero();
  void invert();
  void axpy(const double &, const Departures &);
  double dot_product_with(const Departures &) const;

/// Save departures values
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;

  std::vector<boost::shared_ptr<ObsVector_> > dep_;
};

// =============================================================================

template<typename MODEL>
Departures<MODEL>::Departures(const ObsSpaces_ & obsgeom, const ObsOperators_ & hop): dep_(0)
{
  for (std::size_t jj = 0; jj < obsgeom.size(); ++jj) {
    boost::shared_ptr<ObsVector_> tmp(new ObsVector_(obsgeom[jj], hop[jj].observed()));
    dep_.push_back(tmp);
  }
  Log::trace() << "Departures created" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL>::Departures(std::vector<boost::shared_ptr<ObsVector_> > val)
  : dep_(val)
{
  Log::trace() << "Departures created" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL>::Departures(const Departures & other): dep_(0)
{
// We want deep copy here
  for (std::size_t jj = 0; jj < other.dep_.size(); ++jj) {
    boost::shared_ptr<ObsVector_> tmp(new ObsVector_(*other.dep_[jj]));
    dep_.push_back(tmp);
  }
  Log::trace() << "Departures copy-created" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL>::~Departures() {
  Log::trace() << "Departures destructed" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator=(const Departures & rhs) {
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    *dep_[jj] = *rhs.dep_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator+=(const Departures & rhs) {
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    *dep_[jj] += *rhs.dep_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator-=(const Departures & rhs) {
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    *dep_[jj] -= *rhs.dep_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator*=(const double & zz) {
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    *dep_[jj] *= zz;
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator*=(const Departures & rhs) {
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    *dep_[jj] *= *rhs.dep_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator/=(const Departures & rhs) {
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    *dep_[jj] /= *rhs.dep_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::zero() {
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj]->zero();
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::invert() {
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj]->invert();
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::axpy(const double & zz, const Departures & rhs) {
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj]->axpy(zz, *rhs.dep_[jj]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double Departures<MODEL>::dot_product_with(const Departures & other) const {
  double zz = 0.0;
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    zz += dot_product(*dep_[jj], *other.dep_[jj]);
  }
  return zz;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Departures<MODEL>::save(const std::string & name) const {
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj]->save(name);
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Departures<MODEL>::print(std::ostream & os) const {
  for (std::size_t jj = 0; jj < dep_.size(); ++jj) {
    os << *dep_[jj] << std::endl;
  }
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_DEPARTURES_H_
