/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "atlas/field.h"
#include "oops/util/Printable.h"


namespace oops {

class Variables;

// -----------------------------------------------------------------------------

class LocalInterpolatorBase : public util::Printable {
 public:
  LocalInterpolatorBase() = default;
  virtual ~LocalInterpolatorBase() = default;

  /// Interpolate Variables from source fields to target fields (no
  /// mask).
  virtual void apply(const Variables& variables, const atlas::FieldSet& sourceFieldSet,
                     std::vector<double>& targetFieldVec) const = 0;

  /// Interpolate Variables from source fields to target fields.
  virtual void apply(const Variables& variables, const atlas::FieldSet& sourceFieldSet,
                     const std::vector<bool>& mask,
                     std::vector<double>& targetFieldsVec) const = 0;

  /// Adjoint of interpolation from source to target fields (no mask).
  virtual void applyAD(const Variables& variables, atlas::FieldSet& sourceFieldSet,
                       const std::vector<double>& targetFieldVec) const = 0;

  /// Adjoint of interpolation from source to target fields.
  virtual void applyAD(const Variables& variables, atlas::FieldSet& sourceFieldSet,
                       const std::vector<bool>& mask,
                       const std::vector<double>& targetFieldVec) const = 0;

  // Unscramble MPI buffer into the model's FieldSet representation.
  // Methods do NOT rely on any internal state of the interpolator, they only encode the
  // inverse of the transformation done in apply() to get an MPI buffer from the FieldSet
  static void bufferToFieldSet(const Variables &, const std::vector<size_t> &,
                                const std::vector<double> &, atlas::FieldSet &);
  static void bufferToFieldSetAD(const Variables &, const std::vector<size_t> &,
                                  std::vector<double> &, const atlas::FieldSet &);
};

// -----------------------------------------------------------------------------

}  // namespace oops
