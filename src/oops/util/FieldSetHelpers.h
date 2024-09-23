/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/Variables.h"

#include "eckit/mpi/Comm.h"

namespace util {

enum class ToleranceType {absolute, relative, normalized_absolute};

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

/// Depending on toltype it will compare the values of two fieldsets
/// and output boolean true if the difference is smaller than the tolerance allowed.
/// There are 3 different tolerance conditions:
/// Tolerance::absolute (default) - assumes that the difference is smaller than a
///                    prescribed value (like 1e-4) for all values and all fields.
/// Tolerance::normalized_absolute (similar to Tolerance::absolute but the fixed tolerance
///                                 value is multiplied by the absolute maximum value on
///                                 the PE for each field and model level)
/// Tolerance::relative - each value in the fieldset is compared element-wise.
///                       To be true the absolute difference between two values is the same
///                       or less than the maximum absolute value of the two values multiplied
///                       by the tolerance value.
bool compareFieldSets(const eckit::mpi::Comm & comm,
                      const atlas::FieldSet & fset1,
                      const atlas::FieldSet & fset2,
                      const double & tol = 1.0e-12,
                      const ToleranceType & toltype = ToleranceType::absolute);

bool compareFieldSets(const atlas::FieldSet &,
                      const atlas::FieldSet &,
                      const double & tol = 1.0e-12,
                      const bool & absolute = true);
std::string getGridUid(const atlas::FunctionSpace &);
std::string getGridUid(const atlas::FieldSet &);

/// Extract coordinates and values of Field points in dataFieldSet at indices
/// where the corresponding Field point in diagFieldSet is unity
std::tuple< std::vector<double>,
            std::vector<double>,
            std::vector<size_t>,
            std::vector<size_t>,
            std::vector<double>,
            std::vector<size_t>>
extractUnityPoints(const eckit::mpi::Comm &,
                   const eckit::mpi::Comm &,
                   const atlas::FunctionSpace &,
                   const atlas::FieldSet & dataFieldSet,
                   const atlas::FieldSet & diagFieldSet);

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
