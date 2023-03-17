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

#include "eckit/mpi/Comm.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

namespace util {

// -----------------------------------------------------------------------------

atlas::FieldSet createRandomFieldSet(const oops::GeometryData &,
                                     const std::vector<size_t> &,
                                     const oops::Variables &,
                                     const size_t & timeRank = 0);
atlas::FieldSet copyFieldSet(const atlas::FieldSet &);
atlas::FieldSet shareFields(const atlas::FieldSet &);
void removeFieldsFromFieldSet(atlas::FieldSet &,
                              const oops::Variables &);
std::string getGridUid(const atlas::FunctionSpace &);
std::string getGridUid(const atlas::FieldSet &);
void zeroFieldSet(atlas::FieldSet &);
void addFieldSets(atlas::FieldSet &,
                  const atlas::FieldSet &);
void multiplyFieldSet(atlas::FieldSet &,
                      const double &);
void multiplyFieldSets(atlas::FieldSet &,
                       const atlas::FieldSet &);
double dotProductFieldSets(const atlas::FieldSet &,
                           const atlas::FieldSet &,
                           const oops::Variables &,
                           const eckit::mpi::Comm &);
void divideFieldSets(atlas::FieldSet &,
                     const  atlas::FieldSet &);
void sqrtFieldSet(atlas::FieldSet &);
void printDiagValues(const eckit::mpi::Comm &,
                     const oops::GeometryData &,
                     const atlas::FieldSet &,
                     const atlas::FieldSet &);

// -----------------------------------------------------------------------------

}  // namespace util
