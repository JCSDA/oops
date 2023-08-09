/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/FieldSet3D.h"

#include <cmath>

#include "atlas/array.h"

#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

FieldSet3D::FieldSet3D(const atlas::FieldSet & fset,
                       const util::DateTime & validTime,
                       const eckit::mpi::Comm & comm):
  fset_(fset), validTime_(validTime), vars_(currentVariables()),
  comm_(comm)
{}

// -----------------------------------------------------------------------------

oops::Variables FieldSet3D::currentVariables() const {
  oops::Variables vars;
  for (const auto & field : fset_) {
    vars.push_back(field.name());
    vars.addMetaData(field.name(), "levels", field.levels());
  }
  return vars;
}

// -----------------------------------------------------------------------------

const oops::Variables & FieldSet3D::variables() const {
  vars_ = currentVariables();
  return vars_;
}

// -----------------------------------------------------------------------------

void FieldSet3D::print(std::ostream & os) const {
  os << std::endl << "Valid time: " << validTime() << std::endl;
  for (const auto & field : fset_) {
    os << " Field " << field.name()
       << " norm: " << util::normField(field, comm_) << std::endl;
  }
}

}  // namespace oops
