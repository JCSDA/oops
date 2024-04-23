/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/FieldSets.h"

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

FieldSets::FieldSets(const atlas::FunctionSpace & fspace,
                     const Variables & vars,
                     const std::vector<util::DateTime> & times,
                     const eckit::Configuration & config,
                     const eckit::mpi::Comm & commGeom,
                     const eckit::mpi::Comm & commTime,
                     const eckit::mpi::Comm & commEns):
  Base_(commTime, commEns) {
  Log::trace() << "FieldSets::FieldSets read start " << config << std::endl;

  std::vector<eckit::LocalConfiguration> locals = this->configure(config);

  size_t mytime = this->local_time_size() * commTime.rank();
  size_t indx = 0;
  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      this->dataset().emplace_back(std::make_unique<FieldSet3D>(
                                       times[mytime + jt],
                                       commGeom));
      this->dataset().back()->read(fspace, vars, locals.at(indx));
      ++indx;
    }
  }

  this->sync_times();
  this->check_consistency();

  Log::trace() << "FieldSets::FieldSets read done" << std::endl;
}

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
