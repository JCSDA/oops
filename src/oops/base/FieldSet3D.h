/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <vector>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Wrapper around atlas::FieldSet, storing the valid time associated with the fieldset.
/// FieldSet3D implements all methods that are required to be implemented in the DATA
/// associated with oops::DataSetBase<DATA,...> (base class for the 4D ensemble storage).
/// Note: FieldSet3D assumes that the valid time of the data does not change during
/// the object's lifecycle, unlike oops::State and oops::Increment.
/// FieldSet3D assumes that the variables of the data can change during its lifecycle.
/// Copies of the objects of this class do shallow copies of atlas fieldsets.
class FieldSet3D : public util::Serializable,
                   public util::Printable {
 public:
  FieldSet3D(const atlas::FieldSet & fset, const util::DateTime & validTime,
             const eckit::mpi::Comm & comm);
  /// @brief Creates FieldSet3D with an empty atlas::FieldSet
  FieldSet3D(const util::DateTime & validTime, const eckit::mpi::Comm & comm);

  const Variables & variables() const;
  const util::DateTime validTime() const {return validTime_;}
  const eckit::mpi::Comm & commGeom() const {return comm_;}

  const atlas::FieldSet & fieldSet() const {return fset_;}
  atlas::FieldSet & fieldSet() {return fset_;}

  void zero();
  FieldSet3D & operator+=(const FieldSet3D & other);
  /// @brief Computes dot product of this FieldSet3D with the \p other FieldSet3D
  ///        only for specified variables \p vars.
  double dot_product_with(const FieldSet3D &, const Variables & vars) const;

  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  oops::Variables currentVariables() const;
  void print(std::ostream &) const override;

  atlas::FieldSet fset_;
  const util::DateTime validTime_;
  mutable oops::Variables vars_;
  const eckit::mpi::Comm & comm_;
};

// -----------------------------------------------------------------------------

/// @brief Initializes FieldSet3D to have the same fields as the \p other; values
///        are not copied.
FieldSet3D initFieldSet3D(const FieldSet3D & other);


}  // namespace oops
