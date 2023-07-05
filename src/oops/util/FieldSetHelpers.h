/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/Variables.h"

#include "eckit/mpi/Comm.h"

namespace util {

// -----------------------------------------------------------------------------

atlas::FieldSet createRandomFieldSet(const eckit::mpi::Comm &,
                                     const atlas::FunctionSpace &,
                                     const std::vector<size_t> &,
                                     const std::vector<std::string> &,
                                     const size_t & timeRank = 0);
/// Returns a fieldset with the same smooth field for all variables.
/// Useful for testing interpolation.
atlas::FieldSet createSmoothFieldSet(const eckit::mpi::Comm &,
                                     const atlas::FunctionSpace &,
                                     const std::vector<size_t> &,
                                     const std::vector<std::string> &);

void copyFieldSet(const atlas::FieldSet &, atlas::FieldSet &);
atlas::FieldSet copyFieldSet(const atlas::FieldSet &);
void shareFields(const atlas::FieldSet &, atlas::FieldSet &);
atlas::FieldSet shareFields(const atlas::FieldSet &);
void removeFieldsFromFieldSet(atlas::FieldSet &,
                              const std::vector<std::string> &);
bool compareFieldSets(const atlas::FieldSet &,
                      const atlas::FieldSet &,
                      const double & tol = 1.0e-12);
std::string getGridUid(const atlas::FunctionSpace &);
std::string getGridUid(const atlas::FieldSet &);
void printDiagValues(const eckit::mpi::Comm &,
                     const eckit::mpi::Comm &,
                     const atlas::FunctionSpace &,
                     const atlas::FieldSet &,
                     const atlas::FieldSet &);
void readFieldSet(const eckit::mpi::Comm &,
                  const atlas::FunctionSpace &,
                  const std::vector<size_t> &,
                  const std::vector<std::string> &,
                  const eckit::Configuration &,
                  atlas::FieldSet &);

void writeFieldSet(const eckit::mpi::Comm &,
                   const eckit::Configuration &,
                   const atlas::FieldSet &);

atlas::FieldSet createRandomFieldSet(const eckit::mpi::Comm &,
                                     const atlas::FunctionSpace &,
                                     const oops::Variables &,
                                     const size_t & timeRank = 0);
/// Returns a fieldset with the same smooth field for all variables.
/// Useful for testing interpolation.
atlas::FieldSet createSmoothFieldSet(const eckit::mpi::Comm &,
                                     const atlas::FunctionSpace &,
                                     const oops::Variables &);

void readFieldSet(const eckit::mpi::Comm &,
                  const atlas::FunctionSpace &,
                  const oops::Variables &,
                  const eckit::Configuration &,
                  atlas::FieldSet &);


// -----------------------------------------------------------------------------

}  // namespace util
