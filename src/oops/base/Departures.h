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
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"

namespace oops {

/// Difference between two observation vectors.
/*!
 * A departure is the difference between two observations.
 * The archetypal example is \f$ \mathbf{y} - {\cal H}(\mathbf{x}) \f$.
 *
 * Keeping an observation space vector here is necessary for the implementation
 * of generic observation error covariance matrices.
 */

// -----------------------------------------------------------------------------
template <typename OBS>
class Departures : public GeneralizedDepartures {
  typedef ObsSpaces<OBS>           ObsSpaces_;
  typedef ObsVector<OBS>           ObsVector_;
  template <typename DATA> using ObsData_ = ObsDataVector<OBS, DATA>;
  template <typename DATA> using ObsDataVec_ = std::vector<std::shared_ptr<ObsData_<DATA>>>;

 public:
/// \brief create Departures for all obs (read from ObsSpace if \p name is specified)
  explicit Departures(const ObsSpaces_ &, const std::string & name = "");

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
  void ones();
  void random();
  void invert();
  void axpy(const double &, const Departures &);
  double dot_product_with(const Departures &) const;
  double rms() const;
/// Return number of departures (excluding departures that are masked out)
  size_t nobs() const;

/// Mask out departures where the passed in qc flags are > 0
  void mask(ObsDataVec_<int>);
/// Mask out departures where \p  mask has missing values
  void mask(const Departures & mask);

/// Pack departures in an Eigen vector (excluding departures that are masked out)
  Eigen::VectorXd packEigen(const Departures &) const;
/// Size of departures packed into an Eigen vector
  size_t packEigenSize(const Departures &) const;

/// Save departures values
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;

/// Data
  std::vector<ObsVector_> dep_;
};

// =============================================================================

template<typename OBS>
Departures<OBS>::Departures(const ObsSpaces_ & obsdb,
                            const std::string & name): dep_()
{
  dep_.reserve(obsdb.size());
  for (size_t jj = 0; jj < obsdb.size(); ++jj) {
    dep_.emplace_back(obsdb[jj], name);
  }
  Log::trace() << "Departures created" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename OBS>
Departures<OBS> & Departures<OBS>::operator+=(const Departures & rhs) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj] += rhs[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
Departures<OBS> & Departures<OBS>::operator-=(const Departures & rhs) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj] -= rhs[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
Departures<OBS> & Departures<OBS>::operator*=(const double & zz) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj] *= zz;
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
Departures<OBS> & Departures<OBS>::operator*=(const Departures & rhs) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj] *= rhs[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
Departures<OBS> & Departures<OBS>::operator/=(const Departures & rhs) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj] /= rhs[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void Departures<OBS>::zero() {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj].zero();
  }
}
// -----------------------------------------------------------------------------
template<typename OBS>
void Departures<OBS>::ones() {
  for (auto & dep : dep_) {
    dep.ones();
  }
}
// -----------------------------------------------------------------------------
template<typename OBS>
void Departures<OBS>::random() {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj].random();
  }
}
// -----------------------------------------------------------------------------
template<typename OBS>
void Departures<OBS>::invert() {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj].invert();
  }
}
// -----------------------------------------------------------------------------
template<typename OBS>
void Departures<OBS>::axpy(const double & zz, const Departures & rhs) {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj].axpy(zz, rhs[jj]);
  }
}
// -----------------------------------------------------------------------------
template<typename OBS>
double Departures<OBS>::dot_product_with(const Departures & other) const {
  double zz = 0.0;
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    zz += dot_product(dep_[jj], other[jj]);
  }
  return zz;
}
// -----------------------------------------------------------------------------
template<typename OBS>
double Departures<OBS>::rms() const {
  double zz = 0.0;
  if (nobs() > 0) zz = sqrt(dot_product_with(*this) / this->nobs());
  return zz;
}
// -----------------------------------------------------------------------------
template<typename OBS>
size_t Departures<OBS>::nobs() const {
  size_t nobs = 0;
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    nobs += dep_[jj].nobs();
  }
  return nobs;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void Departures<OBS>::mask(ObsDataVec_<int> qcflags) {
  for (size_t ii = 0; ii < dep_.size(); ++ii) {
    dep_[ii].mask(*qcflags[ii]);
  }
}
// -----------------------------------------------------------------------------
template<typename OBS>
void Departures<OBS>::mask(const Departures & mask) {
  for (size_t ii = 0; ii < dep_.size(); ++ii) {
    dep_[ii].mask(mask[ii]);
  }
}
// -----------------------------------------------------------------------------
template <typename OBS>
Eigen::VectorXd Departures<OBS>::packEigen(const Departures & mask) const {
  std::vector<size_t> len(dep_.size());
  for (size_t idep = 0; idep < dep_.size(); ++idep) {
    len[idep] = dep_[idep].packEigenSize(mask[idep]);
  }
  size_t all_len = std::accumulate(len.begin(), len.end(), 0);

  Eigen::VectorXd vec(all_len);
  size_t ii = 0;
  for (size_t idep = 0; idep < dep_.size(); ++idep) {
    vec.segment(ii, len[idep]) = dep_[idep].packEigen(mask[idep]);
    ii += len[idep];
  }
  return vec;
}
// -----------------------------------------------------------------------------
template <typename OBS>
size_t Departures<OBS>::packEigenSize(const Departures & mask) const {
  size_t len = 0;
  for (size_t idep = 0; idep < dep_.size(); ++idep) {
    len += dep_[idep].packEigenSize(mask[idep]);
  }
  return len;
}
// -----------------------------------------------------------------------------
template <typename OBS>
void Departures<OBS>::save(const std::string & name) const {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    dep_[jj].save(name);
  }
}
// -----------------------------------------------------------------------------
template <typename OBS>
void Departures<OBS>::print(std::ostream & os) const {
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    os << dep_[jj] << std::endl;
  }
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_DEPARTURES_H_
