
// (C) Crown Copyright 2022 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#pragma once

#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/interpolation/Interpolation.h"
#include "eckit/config/Configuration.h"
#include "oops/atlas/Interpolator.h"
#include "oops/base/GeometryData.h"
#include "oops/generic/InterpolationCache.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {

class Variables;

class AtlasInterpolator : public atlasbase::Interpolator,
                          public util::Printable,
                          private util::ObjectCounter<AtlasInterpolator> {
 public:
  /// \brief Class name string.
  static std::string classname() { return "oops::AtlasInterpolator"; }

  /// \brief Construct interpolator from GeometryData and target lat-lons.
  AtlasInterpolator(const eckit::Configuration& conf,
                    const GeometryData& geomData,
                    const std::vector<double>& targetLats,
                    const std::vector<double>& targetLons);

  using atlasbase::Interpolator::apply;
  using atlasbase::Interpolator::applyAD;

  /// \brief Interpolate Variables from source fields to target fields.
  void apply(const Variables& variables, const atlas::FieldSet& sourceFieldSet,
             const std::vector<bool>& mask,
             std::vector<double>& targetFieldsVec) const override;

  /// \brief Adjoint of interpolation from source to target fields.
  void applyAD(const Variables& variables, atlas::FieldSet& sourceFieldSet,
               const std::vector<bool>& mask,
               const std::vector<double>& targetFieldVec) const override;

 private:
  /// \brief Collect handles to source fields in Variables into a new FieldSet.
  atlas::FieldSet copySourceFieldHandles(
      const Variables& variables, const atlas::FieldSet& sourceFields) const;

  /// \brief Create target FieldSet to match source FieldSet.
  atlas::FieldSet createTargetFields(
      const atlas::FunctionSpace& targetFunctionSpace,
      const atlas::FieldSet& sourceFields) const;

  /// \brief Copy interpolation target field set to target field vector.
  void copyData(const atlas::FieldSet& targetFieldSet,
               const std::vector<bool>& mask,
               std::vector<double>& targetFieldVec) const;

  /// \brief Copy target field vector to target field set.
  void copyDataAD(atlas::FieldSet& targetFieldSet,
                 const std::vector<bool>& mask,
                 const std::vector<double>& targetFieldVec) const;

  /// \brief Print object information.
  void print(std::ostream& os) const override;

  // Lazily builds and caches atlas::Interpolation objects, keyed by mask.
  InterpolationCache interpolationCache_;
};

}  // namespace oops
