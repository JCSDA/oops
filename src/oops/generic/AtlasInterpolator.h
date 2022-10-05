
// (C) Crown Copyright 2022 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/util/Point.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredInterpolator.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {

template <typename MODEL>
class AtlasInterpolator
    : public util::Printable,
      private util::ObjectCounter<AtlasInterpolator<MODEL>> {
  typedef Geometry<MODEL> Geometry_;
  typedef State<MODEL> State_;
  typedef Increment<MODEL> Increment_;

 public:
  /// Class name string.
  static std::string classname() { return "oops::AtlasInterpolator"; }

  /// Construct interpolator from Geometry and target lat-lons.
  AtlasInterpolator(const eckit::Configuration& conf, const Geometry_& geom,
                    const std::vector<double>& targetLats,
                    const std::vector<double>& targetLons);

  /// Destructor.
  ~AtlasInterpolator();

  /// Interpolate Variables from State to targetFields (no mask).
  void apply(const Variables& variables, const State_& state,
             std::vector<double>& targetFieldVec) const;

  /// Interpolate Variables from Increment to target fields (no mask).
  void apply(const Variables& variables, const Increment_& increment,
             std::vector<double>& targetFieldVec) const;

  /// Adjoint of interpolation from Increment to target fields (no mask).
  void applyAD(const Variables& variables, Increment_& inc,
               const std::vector<double>& targetFieldVec) const;

  /// Interpolate Variables from State to targetFields.
  void apply(const Variables& variables, const State_& state,
             const std::vector<bool>& mask,
             std::vector<double>& targetFieldVec) const;

  /// Interpolate Variables from Increment to target fields.
  void apply(const Variables& variables, const Increment_& inc,
             const std::vector<bool>& mask,
             std::vector<double>& targetFieldVec) const;

  /// Adjoint of interpolation from Increment to target fields.
  void applyAD(const Variables& variables, Increment_& inc,
               const std::vector<bool>& mask,
               const std::vector<double>& targetFieldVec) const;

  /// Interpolate Variables from source fields to target fields (no
  /// mask).
  void apply(const Variables& variables, const atlas::FieldSet& sourceFieldSet,
             std::vector<double>& targetFieldsVec) const;

  /// Adjoint of interpolation from source to target fields (no mask).
  void applyAD(const Variables& variables, atlas::FieldSet& sourceFieldSet,
               const std::vector<double>& targetFieldVec) const;

  /// Interpolate Variables from source fields to target fields.
  void apply(const Variables& variables, const atlas::FieldSet& sourceFieldSet,
             const std::vector<bool>& mask,
             std::vector<double>& targetFieldsVec) const;

  /// Adjoint of interpolation from source to target fields.
  void applyAD(const Variables& variables, atlas::FieldSet& sourceFieldSet,
               const std::vector<bool>& mask,
               const std::vector<double>& targetFieldVec) const;

  /// Buffer to FieldSet method from UnstructuredInterpolator.
  static constexpr auto bufferToFieldSet =
      &UnstructuredInterpolator<MODEL>::bufferToFieldSet;

  /// Buffer to FieldSet adjoint method from UnstructuredInterpolator.
  static constexpr auto bufferToFieldSetAD =
      &UnstructuredInterpolator<MODEL>::bufferToFieldSetAD;

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

 private:
  void print(std::ostream& os) const;

  // Copy a target FieldSet to a vector.
  template <typename Functor, typename VecIt>
  void fieldSetToVector(const Variables& variables,
                        const std::vector<bool>& mask,
                        atlas::FieldSet& targetFieldSet, VecIt TargetFieldVecIt,
                        const Functor& dataCopy) const;

  // Get the total number of elements to write to the target vector.
  size_t getTotalElements(const Variables& variables,
                          const atlas::FieldSet& inputFields) const;

  // Get or make an interpolation object using targetLonLats and mask.
  const atlas::Interpolation& getInterp(const std::vector<bool>& mask) const;

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

// Recursive ForEach to visit elements of masked vector with
// atlas::ArrayView<Value, Rank>.
// Note: This iterates over array_view(i, j, k) in k varies fastest order.
template <int Rank, int Dim = Rank>
struct ForEach {
  template <typename Value, typename Functor, typename VecIt, typename... Idxs>
  static void apply(const std::vector<bool>& mask,
                    atlas::array::ArrayView<Value, Rank>& targetFieldView,
                    VecIt& targetFieldVecIt, const Functor& dataCopy,
                    Idxs... idxs) {
    // Iterate over dimension Dim of array.
    for (atlas::idx_t idx = 0; idx < targetFieldView.shape(Dim); ++idx) {
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

template <typename MODEL>
AtlasInterpolator<MODEL>::~AtlasInterpolator() {}

template <typename MODEL>
AtlasInterpolator<MODEL>::AtlasInterpolator(
    const eckit::Configuration& conf, const Geometry_& geom,
    const std::vector<double>& targetLats,
    const std::vector<double>& targetLons)
    : sourceFunctionSpace_{geom.functionSpace()},
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

template <typename MODEL>
void AtlasInterpolator<MODEL>::apply(
    const Variables& variables, const State_& state,
    std::vector<double>& targetFieldVec) const {
  apply(variables, state.stateFields(),
        std::vector<bool>(targetLonLats_.size(), true), targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::apply(
    const Variables& variables, const Increment_& increment,
    std::vector<double>& targetFieldVec) const {
  apply(variables, increment.incrementFields(),
        std::vector<bool>(targetLonLats_.size(), true), targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::applyAD(
    const Variables& variables, Increment_& increment,
    const std::vector<double>& targetFieldVec) const {
  applyAD(variables, increment.incrementFields(),
          std::vector<bool>(targetLonLats_.size(), true), targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::apply(
    const Variables& variables, const Increment_& increment,
    const std::vector<bool>& mask, std::vector<double>& targetFieldVec) const {
  apply(variables, increment.incrementFields(), mask, targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::applyAD(
    const Variables& variables, Increment_& increment,
    const std::vector<bool>& mask,
    const std::vector<double>& targetFieldVec) const {
  applyAD(variables, increment.incrementFields(), mask, targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::apply(
    const Variables& variables, const State_& state,
    const std::vector<bool>& mask, std::vector<double>& targetFieldVec) const {
  apply(variables, state.stateFields(), mask, targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::apply(
    const Variables& variables, const atlas::FieldSet& sourceFieldSet,
    std::vector<double>& targetFieldVec) const {
  apply(variables, sourceFieldSet,
        std::vector<bool>(targetLonLats_.size(), true), targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::apply(
    const Variables& variables, const atlas::FieldSet& sourceFieldSet,
    const std::vector<bool>& mask, std::vector<double>& targetFieldVec) const {

  Log::trace() << classname() + "::apply start" << std::endl;
  util::Timer timer(classname(), "apply");

  // Exit if all mask elements are false.
  if (std::none_of(mask.cbegin(), mask.cend(),
                   [](const bool& maskElem)->bool { return maskElem; })) {
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
  for (const auto& variable : interpVars.variables()) {
    interp.execute(tempSourceFieldSet[variable], targetFieldSet[variable]);
  }

  // Post-process fields.
  postProcessFields(targetFieldSet, mask);

  // Resize targetFieldVec (just in case);
  targetFieldVec.resize(getTotalElements(variables, sourceFieldSet));

  // Copy targetFieldSet to vector.
  const auto dataCopy = [](const double& fieldElem, double& vecElem)->void {
    vecElem = fieldElem;
  };
  fieldSetToVector(variables, mask, targetFieldSet, targetFieldVec.begin(),
                   dataCopy);

  Log::trace() << classname() + "::apply done" << std::endl;
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::applyAD(
    const Variables& variables, atlas::FieldSet& sourceFieldSet,
    const std::vector<double>& targetFieldVec) const {
  applyAD(variables, sourceFieldSet,
          std::vector<bool>(targetLonLats_.size(), true), targetFieldVec);
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::applyAD(
    const Variables& variables, atlas::FieldSet& sourceFieldSet,
    const std::vector<bool>& mask,
    const std::vector<double>& targetFieldVec) const {

  Log::trace() << classname() + "::applyAD start" << std::endl;
  util::Timer timer(classname(), "applyAD");

  // Exit if all mask elements are false.
  if (std::none_of(mask.cbegin(), mask.cend(),
                   [](const bool& maskElem)->bool { return maskElem; })) {
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
  const auto dataCopy = [](double& fieldElem, const double& vecElem)->void {
    fieldElem += vecElem;
  };
  fieldSetToVector(variables, mask, targetFieldSet, targetFieldVec.cbegin(),
                   dataCopy);

  // Post-process fields.
  postProcessFieldsAD(targetFieldSet, mask);

  // Interpolation adjoint.
  const auto interpVars = createInterpVariables(variables);
  for (const auto& variable : interpVars.variables()) {
    interp.execute_adjoint(tempSourceFieldSet[variable],
                           targetFieldSet[variable]);
  }

  // Pre-process fields.
  preProcessFieldsAD(tempSourceFieldSet);

  Log::trace() << classname() + "::applyAD done" << std::endl;
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::preProcessFields(atlas::FieldSet& sourceFields)
    const {
  // Do nothing in base class.
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::preProcessFieldsAD(atlas::FieldSet& sourceFields)
    const {
  // Do nothing in base class.
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::postProcessFields(
    atlas::FieldSet& targetFields, const std::vector<bool>& mask) const {
  // Do nothing in base class.
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::postProcessFieldsAD(
    atlas::FieldSet& targetFields, const std::vector<bool>& mask) const {
  // Do nothing in base class.
}

template <typename MODEL>
atlas::FieldSet AtlasInterpolator<MODEL>::copySourceFields(
    const Variables& variables, const atlas::FieldSet& sourceFieldSet) const {
  auto copiedSourceFieldSet = atlas::FieldSet{};

  // Create new field set based on variables.
  for (const auto& variable : variables.variables()) {
    copiedSourceFieldSet.add(sourceFieldSet[variable]);
  }

  return copiedSourceFieldSet;
}

template <typename MODEL>
atlas::FieldSet AtlasInterpolator<MODEL>::createTargetFields(
    const Variables& variables, const atlas::FunctionSpace& targetFunctionSpace,
    const atlas::FieldSet& sourceFieldSet) const {
  auto targetFieldSet = atlas::FieldSet{};

  // Make new target fields which match source fields.
  for (const auto& variable : variables.variables()) {
    // Get source field.
    const auto& sourceField = sourceFieldSet[variable];

    // Configure field using sourceField properties.
    const auto targetConfig = atlas::option::name(sourceField.name()) |
                              atlas::option::levels(sourceField.levels()) |
                              atlas::option::variables(sourceField.variables());

    auto targetField = targetFieldSet.add(
        targetFunctionSpace.createField<double>(targetConfig));
    zeroField(targetField);
  }

  return targetFieldSet;
}

template <typename MODEL>
Variables AtlasInterpolator<MODEL>::createInterpVariables(
    const Variables& inputVariables) const {
  return inputVariables;
}

template <typename MODEL>
void AtlasInterpolator<MODEL>::zeroField(atlas::Field& field) {
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

template <typename MODEL>
const atlas::Interpolation& AtlasInterpolator<MODEL>::getInterp(
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
    atlas::idx_t idx = 0;
    for (const auto& maskElem : mask) {
      if (maskElem) {
        lonLatView(idx, 0) = targetLonLats_[idx].lon();
        lonLatView(idx, 1) = targetLonLats_[idx].lat();
        idx++;
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

template <typename MODEL>
template <typename Functor, typename VecIt>
void AtlasInterpolator<MODEL>::fieldSetToVector(const Variables& variables,
                                                const std::vector<bool>& mask,
                                                atlas::FieldSet& targetFieldSet,
                                                VecIt TargetFieldVecIt,
                                                const Functor& dataCopy) const {
  for (const auto& variable : variables.variables()) {
    auto targetField = targetFieldSet[variable];

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

template <typename MODEL>
void AtlasInterpolator<MODEL>::print(std::ostream& os) const {
  os << classname() << "<" << MODEL::name() << ">";
}

template <typename MODEL>
size_t AtlasInterpolator<MODEL>::getTotalElements(
    const Variables& variables, const atlas::FieldSet& sourceFields) const {
  size_t totalElements = 0;

  // Loop over fields.
  for (const auto& variable : variables.variables()) {
    const auto field = sourceFields[variable];

    size_t elementsPerLocation = 1;
    // Loop over outer elements of field shape (excluding dim 0, the number of
    // source functionspace points).
    for (atlas::idx_t dim = 1; dim < field.rank(); ++dim) {
      elementsPerLocation *= field.shape(dim);
    }
    totalElements += elementsPerLocation;
  }
  return totalElements * targetLonLats_.size();
}

}  // namespace oops
