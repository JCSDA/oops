/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Wrapper around atlas::FieldSet, storing the valid time associated with the fieldset.
/// FieldSet3D implements all methods that are required to be implemented in the DATA
/// associated with oops::DataSetBase<DATA,...> (base class for the 4D ensemble storage).
/// Note: FieldSet3D assumes that the valid time of the data does not change during
/// the object's lifecycle, unlike oops::State and oops::Increment.
/// FieldSet3D assumes that the variables of the data can change during its lifecycle.
class FieldSet3D : public util::Printable {
 public:
  FieldSet3D(const atlas::FieldSet & fset, const util::DateTime & validTime,
             const eckit::mpi::Comm & comm);

  const Variables & variables() const;
  const util::DateTime validTime() const {return validTime_;}
  const atlas::FieldSet & fieldSet() const {return fset_;}
  atlas::FieldSet & fieldSet() {return fset_;}

 private:
  oops::Variables currentVariables() const;
  void print(std::ostream &) const override;

  atlas::FieldSet fset_;
  const util::DateTime validTime_;
  mutable oops::Variables vars_;
  const eckit::mpi::Comm & comm_;
};

// -----------------------------------------------------------------------------

}  // namespace oops
