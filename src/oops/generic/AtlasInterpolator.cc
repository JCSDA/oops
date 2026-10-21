
// (C) Crown Copyright 2022-2026 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#include "oops/generic/AtlasInterpolator.h"

#include <string>

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/option.h"
#include "oops/base/Variables.h"
#include "oops/generic/MaskedVectorFieldSet.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {

AtlasInterpolator::AtlasInterpolator(const eckit::Configuration& conf,
                                     const GeometryData& geomData,
                                     const std::vector<double>& targetLats,
                                     const std::vector<double>& targetLons)
    : atlasbase::Interpolator(conf, geomData, targetLats, targetLons),
      interpolationCache_{geomData.functionSpace(),
                          conf.getSubConfiguration("interpolation method"),
                          targetLats, targetLons} {
  Log::trace() << classname() + "::AtlasInterpolator start" << std::endl;
  Log::trace() << classname() + "::AtlasInterpolator done" << std::endl;
}

void AtlasInterpolator::apply(const Variables& variables,
                              const atlas::FieldSet& sourceFieldSet,
                              const std::vector<bool>& mask,
                              std::vector<double>& targetFieldVec) const {
  Log::trace() << classname() + "::apply start" << std::endl;
  util::Timer timer(classname(), "apply");

  // Get atlas interpolation object.
  const auto& interp = interpolationCache_.getInterp(mask);

  // Get reduced set of source fields.
  const auto interpSourceFieldSet =
      copySourceFieldHandles(variables, sourceFieldSet);

  // Create target fields from reduced source fields.
  auto targetFieldSet =
      createTargetFields(interp.target(), interpSourceFieldSet);

  // Perform interpolation.
  interp.execute(interpSourceFieldSet, targetFieldSet);

  // Copy interpolated data to target field vector.
  copyData(targetFieldSet, mask, targetFieldVec);

  Log::trace() << classname() + "::apply done" << std::endl;
}

void AtlasInterpolator::applyAD(
    const Variables& variables, atlas::FieldSet& sourceFieldSet,
    const std::vector<bool>& mask,
    const std::vector<double>& targetFieldVec) const {
  Log::trace() << classname() + "::applyAD start" << std::endl;
  util::Timer timer(classname(), "applyAD");

  // Get atlas interpolation object.
  const auto& interp = interpolationCache_.getInterp(mask);

  // Get reduced set of source fields.
  auto interpSourceFieldSet = copySourceFieldHandles(variables, sourceFieldSet);

  // Create target fields from reduced source fields.
  auto targetFieldSet =
      createTargetFields(interp.target(), interpSourceFieldSet);

  // Copy target field vector to target field set.
  copyDataAD(targetFieldSet, mask, targetFieldVec);

  // Perform interpolation adjoint.
  interp.execute_adjoint(interpSourceFieldSet, targetFieldSet);

  Log::trace() << classname() + "::applyAD done" << std::endl;
}

atlas::FieldSet AtlasInterpolator::copySourceFieldHandles(
    const Variables& variables, const atlas::FieldSet& sourceFieldSet) const {
  Log::trace() << classname() + "::copySourceFieldHandles start" << std::endl;
  util::Timer timer(classname(), "copySourceFieldHandles");

  auto copiedSourceFieldSet = atlas::FieldSet{};

  // Create new field set based on variables.
  for (const auto& variable : variables) {
    copiedSourceFieldSet.add(sourceFieldSet[variable.name()]);
  }

  Log::trace() << classname() + "::copySourceFieldHandles done" << std::endl;
  return copiedSourceFieldSet;
}

atlas::FieldSet AtlasInterpolator::createTargetFields(
    const atlas::FunctionSpace& targetFunctionSpace,
    const atlas::FieldSet& sourceFieldSet) const {
  Log::trace() << classname() + "::createTargetFields start" << std::endl;
  util::Timer timer(classname(), "createTargetFields");

  auto targetFieldSet = atlas::FieldSet{};

  // Make new target fields which match source fields.
  for (const auto& sourceField : sourceFieldSet) {
    auto targetField = targetFieldSet.add(
        targetFunctionSpace.createField<double>(sourceField.metadata()));
  }

  Log::trace() << classname() + "::createTargetFields done" << std::endl;
  return targetFieldSet;
}

void AtlasInterpolator::copyData(const atlas::FieldSet& targetFieldSet,
                                 const std::vector<bool>& mask,
                                 std::vector<double>& targetFieldVec) const {
  Log::trace() << classname() + "::copyData start" << std::endl;
  util::Timer timer(classname(), "copyData");

  fieldSetToVector(mask, targetFieldSet, targetFieldVec);

  Log::trace() << classname() + "::copyData done" << std::endl;
}

void AtlasInterpolator::copyDataAD(
    atlas::FieldSet& targetFieldSet, const std::vector<bool>& mask,
    const std::vector<double>& targetFieldVec) const {
  Log::trace() << classname() + "::copyDataAD start" << std::endl;
  util::Timer timer(classname(), "copyDataAD");

  vectorToFieldSet(mask, targetFieldVec, targetFieldSet);

  Log::trace() << classname() + "::copyDataAD done" << std::endl;
}

void AtlasInterpolator::print(std::ostream& os) const { os << classname(); }

}  // namespace oops
