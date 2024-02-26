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

#include "eckit/mpi/Comm.h"

namespace util {

// -----------------------------------------------------------------------------

void zeroFieldSet(atlas::FieldSet &);
void addFieldSets(atlas::FieldSet &,
                  const atlas::FieldSet &);
void subtractFieldSets(atlas::FieldSet &,
                       const atlas::FieldSet &);
void multiplyFieldSet(atlas::FieldSet &,
                      const double &);
void multiplyFieldSets(atlas::FieldSet &,
                       const atlas::FieldSet &);
double dotProductFields(const atlas::Field &,
                        const atlas::Field &,
                        const eckit::mpi::Comm &,
                        const bool & includeHalo = true);
double dotProductFieldSets(const atlas::FieldSet &,
                           const atlas::FieldSet &,
                           const std::vector<std::string> &,
                           const eckit::mpi::Comm &,
                           const bool & includeHalo = true);
double normField(const atlas::Field &,
                 const eckit::mpi::Comm &);
double normFieldSet(const atlas::FieldSet &,
                    const std::vector<std::string> &,
                    const eckit::mpi::Comm &);
void divideFieldSets(atlas::FieldSet &,
                     const  atlas::FieldSet &);
void divideFieldSets(atlas::FieldSet &,
                     const  atlas::FieldSet &,
                     const  atlas::FieldSet &);
void sqrtFieldSet(atlas::FieldSet &);
void addZeroFieldToFieldSet(const std::string &,
                            const std::string &,
                            atlas::FieldSet & fset);


// -----------------------------------------------------------------------------

}  // namespace util
