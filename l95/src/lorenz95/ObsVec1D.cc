/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ObsVec1D.h"

#include <math.h>

#include <boost/foreach.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/ObsTableView.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace lorenz95 {
// -----------------------------------------------------------------------------
ObsVec1D::ObsVec1D(const ObsTableView & ot,
                   const std::string & name, const bool fail)
  : obsdb_(ot), data_(ot.nobs())
{
  BOOST_FOREACH(double & val, data_) val = 0.0;
  if (!name.empty()) {
    if (fail || obsdb_.has(name)) obsdb_.getdb(name, data_);
  }
}
// -----------------------------------------------------------------------------
ObsVec1D::ObsVec1D(const ObsVec1D & other)
  : obsdb_(other.obsdb_), data_(other.data_.size())
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
  BOOST_FOREACH(double & val, data_) val *= zz;
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec1D & ObsVec1D::operator+= (const ObsVec1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  for (unsigned int jj = 0; jj < data_.size(); ++jj) data_[jj] += rhs.data_[jj];
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec1D & ObsVec1D::operator-= (const ObsVec1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  for (unsigned int jj = 0; jj < data_.size(); ++jj) data_[jj] -= rhs.data_[jj];
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec1D & ObsVec1D::operator*= (const ObsVec1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  for (unsigned int jj = 0; jj < data_.size(); ++jj) data_[jj] *= rhs.data_[jj];
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec1D & ObsVec1D::operator/= (const ObsVec1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  for (unsigned int jj = 0; jj < data_.size(); ++jj) data_[jj] /= rhs.data_[jj];
  return *this;
}
// -----------------------------------------------------------------------------
void ObsVec1D::zero() {
  BOOST_FOREACH(double & val, data_) val = 0.0;
}
// -----------------------------------------------------------------------------
void ObsVec1D::invert() {
  BOOST_FOREACH(double & val, data_) val = 1.0/val;
}
// -----------------------------------------------------------------------------
void ObsVec1D::axpy(const double & zz, const ObsVec1D & rhs) {
  ASSERT(data_.size() == rhs.data_.size());
  for (unsigned int jj = 0; jj < data_.size(); ++jj) data_[jj] += zz * rhs.data_[jj];
}
// -----------------------------------------------------------------------------
void ObsVec1D::random() {
  obsdb_.random(data_);
}
// -----------------------------------------------------------------------------
double ObsVec1D::dot_product_with(const ObsVec1D & other) const {
  ASSERT(data_.size() == other.data_.size());
  double zz = 0.0;
  for (unsigned int jj = 0; jj < data_.size(); ++jj) zz += data_[jj] * other.data_[jj];
  return zz;
}
// -----------------------------------------------------------------------------
double ObsVec1D::rms() const {
  double zz = 0.0;
  for (unsigned int jj = 0; jj < data_.size(); ++jj) zz += data_[jj] * data_[jj];
  zz = sqrt(zz/data_.size());
  return zz;
}
// -----------------------------------------------------------------------------
void ObsVec1D::save(const std::string & name) const {
  obsdb_.putdb(name, data_);
}
// -----------------------------------------------------------------------------
void ObsVec1D::print(std::ostream & os) const {
  ASSERT(data_.size() > 0);
  double zmin = data_[0];
  double zmax = data_[0];
  double zavg = 0.0;
  BOOST_FOREACH(const double & val, data_) {
    if (val < zmin) zmin = val;
    if (val > zmax) zmax = val;
    zavg += val;
  }
  zavg /= data_.size();
  os << "Lorenz 95 nobs= " << data_.size() << " Min=" << zmin << ", Max=" << zmax
     << ", Average=" << zavg;
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95
