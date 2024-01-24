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

atlas::FieldSet createFieldSet(const atlas::FunctionSpace &,
                               const std::vector<size_t> &,
                               const std::vector<std::string> &);
atlas::FieldSet createFieldSet(const atlas::FunctionSpace &,
                               const std::vector<size_t> &,
                               const std::vector<std::string> &,
                               const double &);
atlas::FieldSet createRandomFieldSet(const eckit::mpi::Comm &,
                                     const atlas::FunctionSpace &,
                                     const std::vector<size_t> &,
                                     const std::vector<std::string> &);
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
                      const double & tol = 1.0e-12,
                      const bool & absolute = true);
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

void readRank3FieldSet(const atlas::FunctionSpace &,
                       const std::vector<size_t> &,
                       const std::vector<std::string> &,
                       atlas::FieldSet &,
                       const std::string &);

void writeFieldSet(const eckit::mpi::Comm &,
                   const eckit::Configuration &,
                   const atlas::FieldSet &);

void writeRank3FieldSet(const atlas::FieldSet &,
                        const std::vector<std::string> &,
                        const atlas::FunctionSpace &,
                        const std::string &,
                        const double &);

atlas::FieldSet createFieldSet(const atlas::FunctionSpace &,
                               const oops::Variables &);
atlas::FieldSet createFieldSet(const atlas::FunctionSpace &,
                               const oops::Variables &,
                               const double &);
atlas::FieldSet createRandomFieldSet(const eckit::mpi::Comm &,
                                     const atlas::FunctionSpace &,
                                     const oops::Variables &);
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
