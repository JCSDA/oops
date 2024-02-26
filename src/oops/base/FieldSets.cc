/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/FieldSets.h"

#include "atlas/array.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

FieldSets::FieldSets(const std::vector<util::DateTime> & times,
                                       const eckit::mpi::Comm & commTime,
                                       const std::vector<int> & members,
                                       const eckit::mpi::Comm & commEns):
  Base_(times, commTime, members, commEns) {}

// -----------------------------------------------------------------------------

FieldSets & FieldSets::operator*=(const oops::FieldSet3D & other) {
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj] *= other;
  }
  return *this;
}

// -----------------------------------------------------------------------------

FieldSets & FieldSets::operator*=(const double zz) {
  for (size_t jj = 0; jj < this->size(); ++jj) {
    (*this)[jj] *= zz;
  }
  return *this;
}

// -----------------------------------------------------------------------------

}  // namespace oops
