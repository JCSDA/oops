/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ObsVec1D.h"

#include <math.h>
#include <limits>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "lorenz95/ObsData1D.h"
#include "lorenz95/ObsTable.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace lorenz95 {
// -----------------------------------------------------------------------------
ObsVec1D::ObsVec1D(const ObsTable & ot,
                   const std::string & name)
  : obsdb_(ot), data_(ot.nobs()), missing_(util::missingValue<double>())
{
  for (double & val : data_) { val = 0.0; }
  if (!name.empty()) obsdb_.getdb(name, data_);
}
// -----------------------------------------------------------------------------
ObsVec1D::ObsVec1D(const ObsVec1D & other)
  : obsdb_(other.obsdb_), data_(other.data_.size()), missing_(util::missingValue<double>())
{
  data_ = other.data_;
}
// -----------------------------------------------------------------------------
ObsVec1D & ObsVec1D::operator= (const ObsVec1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  data_ = rhs.data_;
  return *this;
}

// -----------------------------------------------------------------------------
ObsVec1D & ObsVec1D::operator*= (const double & zz) {
  for (double & val : data_) {
    if (val != missing_) val *= zz;
  }
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec1D & ObsVec1D::operator+= (const ObsVec1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if (data_[jj] == missing_ || rhs.data_[jj] == missing_) {
      data_[jj] = missing_;
    } else {
      data_[jj] += rhs.data_[jj];
    }
  }
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec1D & ObsVec1D::operator-= (const ObsVec1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if (data_[jj] == missing_ || rhs.data_[jj] == missing_) {
      data_[jj] = missing_;
    } else {
      data_[jj] -= rhs.data_[jj];
    }
  }
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec1D & ObsVec1D::operator*= (const ObsVec1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if (data_[jj] == missing_ || rhs.data_[jj] == missing_) {
      data_[jj] = missing_;
    } else {
      data_[jj] *= rhs.data_[jj];
    }
  }
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec1D & ObsVec1D::operator/= (const ObsVec1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if (data_[jj] == missing_ || rhs.data_[jj] == missing_) {
      data_[jj] = missing_;
    } else {
      data_[jj] /= rhs.data_[jj];
    }
  }
  return *this;
}
// -----------------------------------------------------------------------------
void ObsVec1D::zero() {
  for (double & val : data_) val = 0.0;
}
// -----------------------------------------------------------------------------
void ObsVec1D::ones() {
  for (double & val : data_) val = 1.0;
}
// -----------------------------------------------------------------------------
void ObsVec1D::invert() {
  for (double & val : data_) {
    if (val != missing_) val = 1.0/val;
  }
}
// -----------------------------------------------------------------------------
void ObsVec1D::axpy(const double & zz, const ObsVec1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if (data_[jj] == missing_ || rhs.data_[jj] == missing_) {
      data_[jj] = missing_;
    } else {
      data_[jj] += zz * rhs.data_[jj];
    }
  }
}
// -----------------------------------------------------------------------------
void ObsVec1D::random() {
  obsdb_.random(data_);
}
// -----------------------------------------------------------------------------
double ObsVec1D::dot_product_with(const ObsVec1D & other) const {
  ASSERT(data_.size() == other.data_.size());
  double zz = 0.0;
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if ((data_[jj] != missing_) && (other.data_[jj] != missing_)) {
      zz += data_[jj] * other.data_[jj];
    }
  }
  return zz;
}
// -----------------------------------------------------------------------------
double ObsVec1D::rms() const {
  double zz = 0.0;
  double iobs = 0.0;
  for (const double & val : data_) {
    if (val != missing_) {
      zz += val * val;
      iobs += 1.0;
    }
  }
  if (iobs > 0.0) zz = sqrt(zz/iobs);
  return zz;
}
// -----------------------------------------------------------------------------
unsigned int ObsVec1D::nobs() const {
  return data_.size() - std::count(data_.begin(), data_.end(), missing_);
}
// -----------------------------------------------------------------------------
void ObsVec1D::mask(const ObsVec1D & mask) {
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if (mask[jj] == missing_) data_.at(jj) = missing_;
  }
}
// -----------------------------------------------------------------------------
ObsVec1D & ObsVec1D::operator=(const ObsData1D<float> & rhs) {
  const float fmiss = util::missingValue<float>();
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if (rhs[jj] == fmiss) {
      data_.at(jj) = missing_;
    } else {
      data_.at(jj) = static_cast<double>(rhs[jj]);
    }
  }
  return *this;
}
// -----------------------------------------------------------------------------
void ObsVec1D::save(const std::string & name) const {
  obsdb_.putdb(name, data_);
}
// -----------------------------------------------------------------------------
Eigen::VectorXd ObsVec1D::packEigen(const ObsVec1D & mask) const {
  Eigen::VectorXd vec(packEigenSize(mask));
  size_t ii = 0;
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if ((data_[jj] != missing_) && (mask[jj] != missing_)) {
      vec(ii++) = data_[jj];
    }
  }
  return vec;
}
// -----------------------------------------------------------------------------
size_t ObsVec1D::packEigenSize(const ObsVec1D & mask) const {
  size_t ii = 0;
  for (size_t jj = 0; jj < data_.size(); ++jj) {
    if ((data_[jj] != missing_) && (mask[jj] != missing_)) {
      ii++;
    }
  }
  return ii;
}
// -----------------------------------------------------------------------------
void ObsVec1D::read(const std::string & name) {
  obsdb_.getdb(name, data_);
}
// -----------------------------------------------------------------------------
void ObsVec1D::print(std::ostream & os) const {
  double zmin = std::numeric_limits<double>::max();
  double zmax = std::numeric_limits<double>::lowest();
  double zavg = 0.0;
  size_t iobs = 0;
  for (const double & val : data_) {
    if (val != missing_) {
      if (val < zmin) zmin = val;
      if (val > zmax) zmax = val;
      zavg += val;
      ++iobs;
    }
  }
  if (iobs > 0) {
    zavg /= static_cast<double>(iobs);
    os << "Lorenz 95 nobs= " << iobs << " Min=" << zmin << ", Max=" << zmax
       << ", Average=" << zavg;
  } else {
    os << "Lorenz 95 : No observations";
  }
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95
