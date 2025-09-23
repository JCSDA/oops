/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2025 UCAR.
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
#include <sstream>
#include <string>
#include <vector>

#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

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
  template <typename DATA> using ObsDataVec_ = std::vector<ObsDataVector<OBS, DATA>>;

 public:
/// \brief create Departures for all obs (read from ObsSpace if \p name is specified)
  explicit Departures(const ObsSpaces_ &, const std::string & name = "");

/// Access
  size_t size() const {return dep_.size();}
  ObsVector_ & operator[](const size_t ii) {return dep_.at(ii);}
  const ObsVector_ & operator[](const size_t ii) const {return dep_.at(ii);}
  //  const size_t mc(std::vector<double> & missing_value_count) const;
  std::vector<size_t> mc() const;

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
  size_t serialSize() const;
/// Return a vector of serial size at each obsspace
  std::vector<size_t> serialSizes() const;
/// Return the indices of non-masked out elements in ObsVector at each ObsSpace
  std::vector<std::vector<size_t>> maskAndSerialIndices(const Departures &) const;

/// Mask out departures where the passed in qc flags are > 0
  void mask(ObsDataVec_<int>);
/// Mask out departures where \p  mask has missing values
  void mask(const Departures & mask);

/// Pack departures in an Eigen vector (excluding departures that are masked out)
  Eigen::VectorXd packEigen(const Departures &) const;

/// Save departures values
  void save(const std::string &) const;

  std::string info(const std::string & grep = "") const;
  std::string info(const ObsDataVec_<int> &, const std::string & grep = "") const;

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
std::vector<size_t> Departures<OBS>::mc() const {
  std::cout << " dbg: 0";
  std::vector<size_t> missing_value_count;
  missing_value_count.reserve(dep_.size());
  double missing_ = util::missingValue<double>();
  std::cout << " dbg: 1";

  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    std::vector<double> x = dep_[jj].obsvector().data();
    size_t k=0;
    std::cout << " dbg: 2";

    for (size_t i=0;  i < x.size(); ++i){
      if (x[i] != missing_) ++k;
    }
    missing_value_count.push_back(k);
  }
  return missing_value_count;
}

// -----------------------------------------------------------------------------
template<typename OBS>
size_t Departures<OBS>::serialSize() const {
  size_t serialSize = 0;
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    serialSize += dep_[jj].serialSize();
  }
  return serialSize;
}
// -----------------------------------------------------------------------------
template<typename OBS>
std::vector<size_t> Departures<OBS>::serialSizes() const {
  std::vector<size_t> sizevec(dep_.size());
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    sizevec[jj] = dep_[jj].serialSize();
  }
  return sizevec;
}
// -----------------------------------------------------------------------------
template<typename OBS>
std::vector<std::vector<size_t>> Departures<OBS>::maskAndSerialIndices(const Departures & mask)
const {
  std::vector<std::vector<size_t>> indices;
  indices.reserve(dep_.size());
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    std::vector<size_t> indtmp = dep_[jj].maskAndSerialIndices(mask[jj]);
    indices.emplace_back(indtmp);
  }
  return indices;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void Departures<OBS>::mask(ObsDataVec_<int> qcflags) {
  for (size_t ii = 0; ii < dep_.size(); ++ii) {
    dep_[ii].mask(qcflags[ii]);
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
  std::vector<double> valid_values;
  valid_values.reserve(this->serialSize());
  for (size_t idep = 0; idep < dep_.size(); ++idep) {
    dep_[idep].maskAndSerialize(mask[idep], valid_values);
  }
  const Eigen::VectorXd vec = Eigen::Map<Eigen::VectorXd>(valid_values.data(),
                                                          valid_values.size());
  return vec;
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
    os << std::endl << dep_[jj];
  }
}
// -----------------------------------------------------------------------------
template <typename OBS>
std::string Departures<OBS>::info(const std::string & grep) const {
  std::stringstream ss;
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    ss << dep_[jj].info(grep);
  }
  return ss.str();
}
// -----------------------------------------------------------------------------
template <typename OBS>
std::string Departures<OBS>::info(const ObsDataVec_<int> & qcflags,
                                  const std::string & grep) const {
  std::stringstream ss;
  for (size_t jj = 0; jj < dep_.size(); ++jj) {
    ss << dep_[jj].info(grep, qcflags[jj]);
  }
  return ss.str();
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_DEPARTURES_H_
