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

#include <Eigen/Dense>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

template<typename MODEL> class Observations;

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
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
/// \brief create Departures for all obs (read from ObsSpace if name is specified)
  Departures(const ObsSpaces_ &,
             const std::string & name = "", const bool failIfNameNotFound = true);
/// \brief create local Departures
  Departures(const ObsSpaces_ &, const Departures &);

/// Access
  size_t size() const {return dep_.size();}
  ObsVector_ & operator[](const size_t ii) {return dep_.at(ii);}
  const ObsVector_ & operator[](const size_t ii) const {return dep_.at(ii);}

// Linear algebra operators
  Departures & operator+=(const Departures &);
  Departures & operator-=(const Departures &);
  Departures & operator*=(const double &);
  Departures & operator*=(const Departures &);
  Departures & operator/=(const Departures &);
  void zero();
  void random();
  void invert();
  void axpy(const double &, const Departures &);
  double dot_product_with(const Departures &) const;
  double rms() const;
  size_t nobs() const;

/// Pack operators
  Eigen::MatrixXd  packEigen() const;   // pack departures as a 1D Eigen vector

/// Save departures values
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;

/// Data
  std::vector<ObsVector_> dep_;
};

// =============================================================================

template<typename MODEL>
Departures<MODEL>::Departures(const ObsSpaces_ & obsdb,
                              const std::string & name, const bool fail): dep_()
{
  dep_.reserve(obsdb.size());
  for (size_t jj = 0; jj < obsdb.size(); ++jj) {
    dep_.emplace_back(obsdb[jj], name, fail);
  }
  Log::trace() << "Departures created" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL>::Departures(const ObsSpaces_ & obsdb,
                              const Departures & other): dep_() {
  dep_.reserve(obsdb.size());
  for (size_t jj = 0; jj < other.dep_.size(); ++jj) {
    dep_.emplace_back(obsdb[jj], other[jj]);
  }
  Log::trace() << "Local Departures created" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator+=(const Departures & rhs) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj] += rhs[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator-=(const Departures & rhs) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj] -= rhs[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator*=(const double & zz) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj] *= zz;
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator*=(const Departures & rhs) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj] *= rhs[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
Departures<MODEL> & Departures<MODEL>::operator/=(const Departures & rhs) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj] /= rhs[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::zero() {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj].zero();
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::random() {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj].random();
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::invert() {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj].invert();
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void Departures<MODEL>::axpy(const double & zz, const Departures & rhs) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj].axpy(zz, rhs[jj]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double Departures<MODEL>::dot_product_with(const Departures & other) const {
  double zz = 0.0;
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    zz += dot_product(dep_[jj], other[jj]);
  }
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double Departures<MODEL>::rms() const {
  return sqrt(dot_product_with(*this) / this->nobs());
}
// -----------------------------------------------------------------------------
template<typename MODEL>
size_t Departures<MODEL>::nobs() const {
  size_t nobs = 0;
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    nobs += dep_[jj].nobs();
  }
  return nobs;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
Eigen::MatrixXd Departures<MODEL>::packEigen() const {
  Eigen::MatrixXd data1d(1, this->nobs());
  int i = 0;
  for (size_t idep = 0; idep < dep_.size(); ++idep) {
    for (size_t iob = 0; iob < dep_[idep].nobs(); ++iob) {
      data1d(0, i++) = dep_[idep][iob];
    }
  }
  return data1d;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Departures<MODEL>::save(const std::string & name) const {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj].save(name);
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void Departures<MODEL>::print(std::ostream & os) const {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    os << dep_[jj] << std::endl;
  }
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_DEPARTURES_H_
