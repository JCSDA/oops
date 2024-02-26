
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

#include "oops/base/GeometryData.h"
#include "oops/generic/LocalInterpolatorBase.h"
#include "oops/util/ObjectCounter.h"

namespace oops {

class Variables;

template <typename MODEL>
class Geometry;

template <typename MODEL>
class State;

template <typename MODEL>
class Increment;

class AtlasInterpolator : public LocalInterpolatorBase,
                          private util::ObjectCounter<AtlasInterpolator> {
 public:
  /// Class name string.
  static std::string classname() { return "oops::AtlasInterpolator"; }

  /// Construct interpolator from GeometryData and target lat-lons.
  AtlasInterpolator(const eckit::Configuration& conf,
                    const GeometryData& geomData,
                    const std::vector<double>& targetLats,
                    const std::vector<double>& targetLons);

  /// Construct interpolator from Geometry and target lat-lons.
  template <typename MODEL>
  AtlasInterpolator(const eckit::Configuration& conf,
                    const Geometry<MODEL>& geom,
                    const std::vector<double>& targetLats,
                    const std::vector<double>& targetLons);

  /// Destructor.
  ~AtlasInterpolator();

  /// Interpolate Variables from source fields to target fields (no
  /// mask).
  void apply(const Variables& variables, const atlas::FieldSet& sourceFieldSet,
             std::vector<double>& targetFieldVec) const;

  /// Interpolate Variables from source fields to target fields.
  void apply(const Variables& variables, const atlas::FieldSet& sourceFieldSet,
             const std::vector<bool>& mask,
             std::vector<double>& targetFieldsVec) const;

  /// Adjoint of interpolation from source to target fields (no mask).
  void applyAD(const Variables& variables, atlas::FieldSet& sourceFieldSet,
               const std::vector<double>& targetFieldVec) const;

  /// Adjoint of interpolation from source to target fields.
  void applyAD(const Variables& variables, atlas::FieldSet& sourceFieldSet,
               const std::vector<bool>& mask,
               const std::vector<double>& targetFieldVec) const;

  /// Interpolate Variables from State to targetFields (no mask).
  template <typename MODEL>
  void apply(const Variables& variables, const State<MODEL>& state,
             std::vector<double>& targetFieldVec) const;

  /// Interpolate Variables from State to targetFields.
  template <typename MODEL>
  void apply(const Variables& variables, const State<MODEL>& state,
             const std::vector<bool>& mask,
             std::vector<double>& targetFieldVec) const;

  /// Interpolate Variables from Increment to target fields (no mask).
  template <typename MODEL>
  void apply(const Variables& variables, const Increment<MODEL>& increment,
             std::vector<double>& targetFieldVec) const;

  /// Interpolate Variables from Increment to target fields.
  template <typename MODEL>
  void apply(const Variables& variables, const Increment<MODEL>& inc,
             const std::vector<bool>& mask,
             std::vector<double>& targetFieldVec) const;

  /// Adjoint of interpolation from Increment to target fields (no mask).
  template <typename MODEL>
  void applyAD(const Variables& variables, Increment<MODEL>& inc,
               const std::vector<double>& targetFieldVec) const;

  /// Adjoint of interpolation from Increment to target fields.
  template <typename MODEL>
  void applyAD(const Variables& variables, Increment<MODEL>& inc,
               const std::vector<bool>& mask,
               const std::vector<double>& targetFieldVec) const;

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

template <typename MODEL>
AtlasInterpolator::AtlasInterpolator(const eckit::Configuration& conf,
                                     const Geometry<MODEL>& geom,
                                     const std::vector<double>& targetLats,
                                     const std::vector<double>& targetLons)
    : AtlasInterpolator(conf, geom.generic(), targetLats, targetLons) {}

template <typename MODEL>
void AtlasInterpolator::apply(const Variables& variables,
                              const State<MODEL>& state,
                              std::vector<double>& targetFieldVec) const {
  apply(variables, state.fieldSet().fieldSet(),
        std::vector<bool>(targetLonLats_.size(), true), targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator::apply(const Variables& variables,
                              const State<MODEL>& state,
                              const std::vector<bool>& mask,
                              std::vector<double>& targetFieldVec) const {
  apply(variables, state.fieldSet().fieldSet(), mask, targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator::apply(const Variables& variables,
                              const Increment<MODEL>& increment,
                              const std::vector<bool>& mask,
                              std::vector<double>& targetFieldVec) const {
  apply(variables, increment.fieldSet().fieldSet(), mask, targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator::apply(const Variables& variables,
                              const Increment<MODEL>& increment,
                              std::vector<double>& targetFieldVec) const {
  apply(variables, increment.fieldSet().fieldSet(),
        std::vector<bool>(targetLonLats_.size(), true), targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator::applyAD(
    const Variables& variables, Increment<MODEL>& increment,
    const std::vector<double>& targetFieldVec) const {
  applyAD(variables, increment.fieldSet().fieldSet(),
          std::vector<bool>(targetLonLats_.size(), true), targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator::applyAD(
    const Variables& variables, Increment<MODEL>& increment,
    const std::vector<bool>& mask,
    const std::vector<double>& targetFieldVec) const {
  applyAD(variables, increment.fieldSet().fieldSet(), mask, targetFieldVec);
}
}  // namespace oops
