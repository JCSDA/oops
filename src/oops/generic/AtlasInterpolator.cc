
// (C) Crown Copyright 2022 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#include <iostream>
#include <string>

#include "atlas/util/Point.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/generic/AtlasInterpolator.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {

// Helper structs and functions
namespace {

// Recursive ForEach to visit elements of masked vector with
// atlas::ArrayView<Value, Rank>.
// Note: This iterates over array_view(i, j, k) in i varies fastest order.
template <int Rank, int Dim = Rank>
struct ForEach {
  template <typename Value, typename Functor, typename VecIt, typename... Idxs>
  static void apply(const std::vector<bool>& mask,
                    atlas::array::ArrayView<Value, Rank>& targetFieldView,
                    VecIt& targetFieldVecIt, const Functor& dataCopy,
                    Idxs... idxs) {
    // Iterate over dimension Dim of array.
    for (atlas::idx_t idx = 0; idx < targetFieldView.shape(Dim - 1); ++idx) {
      ForEach<Rank, Dim - 1>::apply(mask, targetFieldView, targetFieldVecIt,
                                    dataCopy, idx, idxs...);
    }
  }
};

// End recursion when Dim == 1
template <int Rank>
struct ForEach<Rank, 1> {
  template <typename Value, typename Functor, typename VecIt, typename... Idxs>
  static void apply(const std::vector<bool>& mask,
                    atlas::array::ArrayView<Value, Rank>& targetFieldView,
                    VecIt& targetFieldVecIt, const Functor& dataCopy,
                    Idxs... idxs) {
    // Iterate over mask and call functor.
    atlas::idx_t idx = 0;
    for (const auto& maskElem : mask) {
      if (maskElem) {
        dataCopy(targetFieldView(idx++, idxs...), *targetFieldVecIt);
      }
      ++targetFieldVecIt;
    }
  }
};

template <typename Functor, typename VecIt>
void fieldSetToVector(const Variables& variables, const std::vector<bool>& mask,
                      atlas::FieldSet& targetFieldSet, VecIt TargetFieldVecIt,
                      const Functor& dataCopy) {
  for (const auto& variable : variables) {
    auto targetField = targetFieldSet[variable.name()];

    switch (targetField.rank()) {
      case 1: {
        auto targetFieldView = atlas::array::make_view<double, 1>(targetField);
        ForEach<1>::apply(mask, targetFieldView, TargetFieldVecIt, dataCopy);
        break;
      }
      case 2: {
        auto targetFieldView = atlas::array::make_view<double, 2>(targetField);
        ForEach<2>::apply(mask, targetFieldView, TargetFieldVecIt, dataCopy);
        break;
      }
      case 3: {
        auto targetFieldView = atlas::array::make_view<double, 3>(targetField);
        ForEach<3>::apply(mask, targetFieldView, TargetFieldVecIt, dataCopy);
        break;
      }
      default: {
        const auto errMsg = "No implementation for rank " +
                            std::to_string(targetField.rank()) + " fields.";
        eckit::NotImplemented(errMsg, Here());
      }
    }
  }
}
}  // namespace

AtlasInterpolator::AtlasInterpolator(const eckit::Configuration& conf,
                                     const GeometryData& geomData,
                                     const std::vector<double>& targetLats,
                                     const std::vector<double>& targetLons)
    : sourceFunctionSpace_{geomData.functionSpace()},
      interpMethod_{conf.getSubConfiguration("interpolation method")} {

  Log::trace() << classname() + "::AtlasInterpolator start" << std::endl;
  util::Timer timer(classname(), "AtlasInterpolator");

  // Normalise and save lon lats.
  targetLonLats_.resize(targetLats.size());
  for (size_t idx = 0; idx < targetLonLats_.size(); ++idx) {
    targetLonLats_[idx] = atlas::PointLonLat(targetLons[idx], targetLats[idx]);
    targetLonLats_[idx].normalise();
  }

  Log::trace() << classname() + "::AtlasInterpolator done" << std::endl;
}

AtlasInterpolator::~AtlasInterpolator() {}

void AtlasInterpolator::apply(const Variables& variables,
                              const atlas::FieldSet& sourceFieldSet,
                              std::vector<double>& targetFieldVec) const {
  apply(variables, sourceFieldSet,
        std::vector<bool>(targetLonLats_.size(), true), targetFieldVec);
}

void AtlasInterpolator::apply(const Variables& variables,
                              const atlas::FieldSet& sourceFieldSet,
                              const std::vector<bool>& mask,
                              std::vector<double>& targetFieldVec) const {
  Log::trace() << classname() + "::apply start" << std::endl;
  util::Timer timer(classname(), "apply");

  // Resize targetFieldVec (just in case);
  targetFieldVec.resize(getTotalElements(variables, sourceFieldSet));

  // Exit if all mask elements are false.
  if (std::none_of(mask.cbegin(), mask.cend(),
                   [](const bool & maskElem)->bool { return maskElem; })) {
    return;
  }

  // Get atlas interpolation object.
  const auto& interp = getInterp(mask);

  // Get reduced set of source fields.
  auto tempSourceFieldSet = copySourceFields(variables, sourceFieldSet);

  // Create target fields from reduced source fields.
  auto targetFieldSet =
      createTargetFields(variables, interp.target(), tempSourceFieldSet);

  // Pre-process fields.
  preProcessFields(tempSourceFieldSet);

  // Perform interpolation.
  const auto interpVars = createInterpVariables(variables);
  for (const auto& variable : interpVars) {
    interp.execute(tempSourceFieldSet[variable.name()], targetFieldSet[variable.name()]);
  }

  // Post-process fields.
  postProcessFields(targetFieldSet, mask);

  // Copy targetFieldSet to vector.
  const auto dataCopy = [](const double & fieldElem, double & vecElem)->void {
    vecElem = fieldElem;
  };
  fieldSetToVector(variables, mask, targetFieldSet, targetFieldVec.begin(),
                   dataCopy);

  Log::trace() << classname() + "::apply done" << std::endl;
}

void AtlasInterpolator::applyAD(
    const Variables& variables, atlas::FieldSet& sourceFieldSet,
    const std::vector<double>& targetFieldVec) const {
  applyAD(variables, sourceFieldSet,
          std::vector<bool>(targetLonLats_.size(), true), targetFieldVec);
}

void AtlasInterpolator::applyAD(
    const Variables& variables, atlas::FieldSet& sourceFieldSet,
    const std::vector<bool>& mask,
    const std::vector<double>& targetFieldVec) const {

  Log::trace() << classname() + "::applyAD start" << std::endl;
  util::Timer timer(classname(), "applyAD");

  // Exit if all mask elements are false.
  if (std::none_of(mask.cbegin(), mask.cend(),
                   [](const bool & maskElem)->bool { return maskElem; })) {
    return;
  }

  // Get atlas interpolation object.
  const auto& interp = getInterp(mask);

  // Get reduced set of source fields.
  auto tempSourceFieldSet = copySourceFields(variables, sourceFieldSet);

  // Create target fields from reduced source fields.
  auto targetFieldSet =
      createTargetFields(variables, interp.target(), tempSourceFieldSet);

  // Copy vector to targetFieldSet.
  const auto dataCopy = [](double & fieldElem, const double & vecElem)->void {
    fieldElem += vecElem;
  };
  fieldSetToVector(variables, mask, targetFieldSet, targetFieldVec.cbegin(),
                   dataCopy);

  // Post-process fields.
  postProcessFieldsAD(targetFieldSet, mask);

  // Interpolation adjoint.
  const auto interpVars = createInterpVariables(variables);
  for (const auto& variable : interpVars) {
    interp.execute_adjoint(tempSourceFieldSet[variable.name()],
                           targetFieldSet[variable.name()]);
  }

  // Pre-process fields.
  preProcessFieldsAD(tempSourceFieldSet);

  Log::trace() << classname() + "::applyAD done" << std::endl;
}

void AtlasInterpolator::preProcessFields(atlas::FieldSet& sourceFields) const {
  // Do nothing in base class.
}

void AtlasInterpolator::preProcessFieldsAD(atlas::FieldSet& sourceFields)
    const {
  // Do nothing in base class.
}

void AtlasInterpolator::postProcessFields(atlas::FieldSet& targetFields,
                                          const std::vector<bool>& mask) const {
  // Do nothing in base class.
}

void AtlasInterpolator::postProcessFieldsAD(
    atlas::FieldSet& targetFields, const std::vector<bool>& mask) const {
  // Do nothing in base class.
}

atlas::FieldSet AtlasInterpolator::copySourceFields(
    const Variables& variables, const atlas::FieldSet& sourceFieldSet) const {
  auto copiedSourceFieldSet = atlas::FieldSet{};

  // Create new field set based on variables.
  for (const auto& variable : variables) {
    copiedSourceFieldSet.add(sourceFieldSet[variable.name()]);
  }
  return copiedSourceFieldSet;
}

atlas::FieldSet AtlasInterpolator::createTargetFields(
    const Variables& variables, const atlas::FunctionSpace& targetFunctionSpace,
    const atlas::FieldSet& sourceFieldSet) const {
  auto targetFieldSet = atlas::FieldSet{};

  // Make new target fields which match source fields.
  for (const auto& variable : variables) {
    // Get source field.
    const auto& sourceField = sourceFieldSet[variable.name()];

    // Configure field using sourceField properties.
    const auto targetConfig = atlas::option::name(sourceField.name()) |
                              atlas::option::levels(sourceField.shape(1)) |
                              atlas::option::variables(sourceField.variables());

    auto targetField = targetFieldSet.add(
        targetFunctionSpace.createField<double>(targetConfig));
    zeroField(targetField);
  }

  return targetFieldSet;
}

Variables AtlasInterpolator::createInterpVariables(
    const Variables& inputVariables) const {
  return inputVariables;
}

void AtlasInterpolator::zeroField(atlas::Field& field) {
  switch (field.rank()) {
    case 1: {
      atlas::array::make_view<double, 1>(field).assign(0.);
      break;
    }
    case 2: {
      atlas::array::make_view<double, 2>(field).assign(0.);
      break;
    }
    case 3: {
      atlas::array::make_view<double, 3>(field).assign(0.);
      break;
    }
    default: {
      const auto errMsg = "No implementation for rank " +
                          std::to_string(field.rank()) + " fields.";
      eckit::NotImplemented(errMsg, Here());
    }
  }
}

const atlas::Interpolation& AtlasInterpolator::getInterp(
    const std::vector<bool>& mask) const {
  // Find or insert interpolation object in map.
  auto& interp = interpMap_[mask];

  if (!interp) {
    // Make a lon-lat field.
    const atlas::idx_t fieldSize = std::count(mask.cbegin(), mask.cend(), true);
    auto lonLatField =
        atlas::Field("lonlat", atlas::array::make_datatype<double>(),
                     atlas::array::make_shape(fieldSize, 2));

    auto lonLatView = atlas::array::make_view<double, 2>(lonLatField);

    // Copy lonlats with mask == true.
    atlas::idx_t maskedIdx = 0;
    for (size_t unmaskedIdx = 0; unmaskedIdx < mask.size(); unmaskedIdx++) {
      if (mask[unmaskedIdx] == true) {
        lonLatView(maskedIdx, 0) = targetLonLats_[unmaskedIdx].lon();
        lonLatView(maskedIdx, 1) = targetLonLats_[unmaskedIdx].lat();
        maskedIdx++;
      }
    }

    // Create an interpolation object.
    const auto targetFunctionSpace =
        atlas::functionspace::PointCloud(lonLatField);
    interp = atlas::Interpolation(interpMethod_, sourceFunctionSpace_,
                                  targetFunctionSpace);
  }
  return interp;
}

void AtlasInterpolator::print(std::ostream& os) const { os << classname(); }

size_t AtlasInterpolator::getTotalElements(
    const Variables& variables, const atlas::FieldSet& sourceFields) const {
  size_t totalElements = 0;

  // Loop over fields.
  for (const auto& variable : variables) {
    const auto field = sourceFields[variable.name()];

    size_t elementsPerLocation = 1;
    // Loop over outer elements of field shape (excluding dim 0, the number of
    // source functionspace points).
    for (atlas::idx_t dim = 1; dim < field.rank(); ++dim) {
      elementsPerLocation *= field.shape(dim);
    }
    totalElements += elementsPerLocation;
  }
  totalElements *= targetLonLats_.size();

  return totalElements;
}

}  // namespace oops
