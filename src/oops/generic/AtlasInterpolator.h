
// (C) Crown Copyright 2022 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/interpolation/Interpolation.h"

#include "eckit/config/Configuration.h"

#include "oops/atlas/Interpolator.h"
#include "oops/base/GeometryData.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {

class Variables;

class AtlasInterpolator : public atlasbase::Interpolator,
                          public util::Printable,
                          private util::ObjectCounter<AtlasInterpolator> {
 public:
  /// Class name string.
  static std::string classname() { return "oops::AtlasInterpolator"; }

  /// Construct interpolator from GeometryData and target lat-lons.
  AtlasInterpolator(const eckit::Configuration& conf,
                    const GeometryData& geomData,
                    const std::vector<double>& targetLats,
                    const std::vector<double>& targetLons);

  /// Destructor.
  ~AtlasInterpolator();

  using atlasbase::Interpolator::apply;
  using atlasbase::Interpolator::applyAD;

  /// Interpolate Variables from source fields to target fields.
  void apply(const Variables& variables, const atlas::FieldSet& sourceFieldSet,
             const std::vector<bool>& mask,
             std::vector<double>& targetFieldsVec) const override;

  /// Adjoint of interpolation from source to target fields.
  void applyAD(const Variables& variables, atlas::FieldSet& sourceFieldSet,
               const std::vector<bool>& mask,
               const std::vector<double>& targetFieldVec) const override;

 protected:
  /// Apply pre-processing to sourceFields (overridable)
  virtual void preProcessFields(atlas::FieldSet& sourceFields) const;

  /// Apply pre-processing adjoint to sourceFields (overridable)
  virtual void preProcessFieldsAD(atlas::FieldSet& sourceFields) const;

  /// Apply post-processing to targetFields (overridable)
  virtual void postProcessFields(atlas::FieldSet& targetFields,
                                 const std::vector<bool>& mask) const;

  /// Apply post-processing adjoint to targetFields (overridable)
  virtual void postProcessFieldsAD(atlas::FieldSet& targetFields,
                                   const std::vector<bool>& mask) const;

  /// Copy source fields in Variables to new FieldSet (overridable).
  virtual atlas::FieldSet copySourceFields(
      const Variables& variables, const atlas::FieldSet& sourceFields) const;

  /// Create target FieldSet to match source FieldSet (overridable).
  virtual atlas::FieldSet createTargetFields(
      const Variables& variables,
      const atlas::FunctionSpace& targetFunctionSpace,
      const atlas::FieldSet& sourceFields) const;

  /// Create Variables object for interpolation call (overridable).
  virtual Variables createInterpVariables(const Variables& inputVariables)
      const;

  /// Set all values in a field to zero.
  static void zeroField(atlas::Field& field);

  /// Get or make an interpolation object using targetLonLats and mask.
  const atlas::Interpolation& getInterp(const std::vector<bool>& mask) const;

 private:
  virtual void print(std::ostream& os) const;

  // Get the total number of elements to write to the target vector.
  size_t getTotalElements(const Variables& variables,
                          const atlas::FieldSet& inputFields) const;

  // FunctionSpace from Geometry.
  atlas::FunctionSpace sourceFunctionSpace_{};

  // map of interpolation objects.
  // Mutable keyword required as we need to construct an interpolation object if
  // it doesn't already exist.
  mutable std::unordered_map<std::vector<bool>, atlas::Interpolation>
      interpMap_{};

  // Vector of lon lats from constructor.
  std::vector<atlas::PointLonLat> targetLonLats_{};

  // Atlas Interpolation method config.
  eckit::LocalConfiguration interpMethod_;
};

}  // namespace oops
