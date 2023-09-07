/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/FieldSet4D.h"

#include <algorithm>
#include <cmath>
#include <memory>

#include "atlas/array.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

FieldSet4D copyFieldSet4D(const FieldSet4D & other) {
  FieldSet4D copy(other.times(), other.commTime(), other[0].commGeom());
  for (size_t jtime = 0; jtime < other.size(); ++jtime) {
    util::copyFieldSet(other[jtime].fieldSet(), copy[jtime].fieldSet());
  }
  return copy;
}

// -----------------------------------------------------------------------------

FieldSet4D::FieldSet4D(const std::vector<util::DateTime> & times,
                       const eckit::mpi::Comm & commTime,
                       const eckit::mpi::Comm & commGeom)
  : DataSetBase<FieldSet3D, atlas::FunctionSpace>(times, commTime, {0}, oops::mpi::myself())
{
  size_t mytime = this->local_time_size() * commTime.rank();
  for (size_t jj = 0; jj < this->local_time_size(); ++jj) {
    this->dataset().push_back(std::make_unique<FieldSet3D>(times[mytime + jj], commGeom));
  }
}

// -----------------------------------------------------------------------------

FieldSet4D::FieldSet4D(const FieldSet3D & fset3d)
  : DataSetBase<FieldSet3D, atlas::FunctionSpace>({fset3d.validTime()}, oops::mpi::myself(),
                                                  {0}, oops::mpi::myself())
{
  this->dataset().push_back(std::make_unique<FieldSet3D>(fset3d));
}

// -----------------------------------------------------------------------------

void FieldSet4D::zero() {
  this->check_consistency();
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj].zero();
  }
}

// -----------------------------------------------------------------------------

FieldSet4D & FieldSet4D::operator+=(const FieldSet4D & other) {
  this->check_consistency(other, false);
  ASSERT(this->is_4d());
  for (size_t jt = 0; jt < this->size(); ++jt) {
    (*this)[jt] += other[jt];
  }
  return *this;
}

// -----------------------------------------------------------------------------

FieldSet4D & FieldSet4D::operator*=(const FieldSet4D & other) {
  this->check_consistency(other, false);
  ASSERT(this->is_4d());
  for (size_t jt = 0; jt < this->size(); ++jt) {
    util::multiplyFieldSets((*this)[jt].fieldSet(), other[jt].fieldSet());
  }
  return *this;
}

// -----------------------------------------------------------------------------

FieldSet4D & FieldSet4D::operator*=(const atlas::FieldSet & other) {
  ASSERT(this->is_4d());
  for (size_t jt = 0; jt < this->size(); ++jt) {
    util::multiplyFieldSets((*this)[jt].fieldSet(), other);
  }
  return *this;
}

// -----------------------------------------------------------------------------

FieldSet4D & FieldSet4D::operator*=(const double zz) {
  ASSERT(this->is_4d());
  for (size_t jj = 0; jj < this->size(); ++jj) {
    util::multiplyFieldSet((*this)[jj].fieldSet(), zz);
  }
  return *this;
}

// -----------------------------------------------------------------------------

double FieldSet4D::dot_product_with(const FieldSet4D & other, const oops::Variables & vars) const {
  this->check_consistency(other);
  ASSERT(this->is_4d());
  double zz = 0.0;
  for (size_t jj = 0; jj < this->size(); ++jj) {
    zz += (*this)[jj].dot_product_with(other[jj], vars);
  }
  this->commTime().allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  return zz;
}

// -----------------------------------------------------------------------------

}  // namespace oops
