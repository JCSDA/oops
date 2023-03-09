
// (C) Crown Copyright 2023 Met Office.

// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#include "oops/generic/HtlmRegularization.h"

#include <numeric>

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "oops/util/abor1_cpp.h"

namespace oops {

HtlmRegularizationPart::HtlmRegularizationPart(const HtlmRegularizationPartParameters & params,
                                               const std::vector<std::string> & allVariables,
                                               const std::vector<size_t> & allLevels,
                                               const std::vector<double> & limitsLon,
                                               const std::vector<double> & limitsLat)
: params_(params), value_(params_.value.value()),
  variables_((params_.variables.value() != boost::none) ?
    *params_.variables.value() : allVariables), levels_(allLevels), boundingLons_(limitsLon),
  boundingLats_(limitsLat), containsAllGridPoints_(true) {
  if (!AIsSubsetOfB(variables_, allVariables)) {
    ABORT("HtlmRegularizationPart: \"variables\" must be a subset of H-TLM training variables");
  }
  if (params_.levels.value() != boost::none) {
    levels_ = *params.levels.value();
    if (!AIsSubsetOfB(levels_, allLevels)) {
      ABORT("HtlmRegularizationPart: \"levels\" must be a subset of model levels indexed from 0");
    }
    containsAllGridPoints_ = false;
  }
  if (params_.boundingLons.value() != boost::none) {
    boundingLons_ = *params.boundingLons.value();
    if (!AllOfAAreInRangeOfB(boundingLons_, limitsLon)) {
      ABORT("HtlmRegularizationPart: \"bounding lons\" must be between -180 and 180 degrees");
    }
    containsAllGridPoints_ = false;
  }
  if (params_.boundingLats.value() != boost::none) {
    boundingLats_ = *params.boundingLats.value();
    if (!AllOfAAreInRangeOfB(boundingLats_, limitsLat)) {
      ABORT("HtlmRegularizationPart: \"bounding lats\" must be between -90 and 90 degrees");
    }
    containsAllGridPoints_ = false;
  }
  if (boundingLons_.size() != 2 || boundingLats_.size() != 2) {
    ABORT("HtlmRegularizationPart: \"bounding lats\"/\"...lons\" must list only 2 values");
  }
}

//--------------------------------------------------------------------------------------------------

const bool HtlmRegularizationPart::AllOfAAreInRangeOfB(const std::vector<double> & A,
                                                       const std::vector<double> & B) const {
  return std::none_of(A.begin(), A.end(),
    [B](double elementOfA){return((*std::min_element(B.begin(), B.end()) > elementOfA)
                                  || (elementOfA > *std::max_element(B.begin(), B.end())));});
}

//--------------------------------------------------------------------------------------------------

HtlmRegularizationComponentDependent::HtlmRegularizationComponentDependent(
                                                        const HtlmRegularizationParameters & params,
                                                        atlas::FieldSet templateFieldSet)
: HtlmRegularization(params), nLevels_(templateFieldSet[0].shape(1)),
  nLocations_(templateFieldSet[0].shape(0)) {
  for (auto & templateField : templateFieldSet) {
    auto templateFieldArray = atlas::array::make_view<double, 2>(templateField);
    templateFieldArray.assign(baseValue_);
  }
  // Generate default HtlmRegularizationPart::(member variables)...
  // ...for ::variables_; all of the variables in templateFieldSet
  std::vector<std::string> allVariables = templateFieldSet.field_names();
  // ...for ::levels_; all of the model levels indexed from 0
  std::vector<size_t> allLevels(nLevels_);
  std::iota(std::begin(allLevels), std::end(allLevels), 0);
  // ...for ::boundingLons_; -180 to 180 degrees
  std::vector<double> limitsLon = {-180.0, 180.0};
  // ...for ::boundingLats_; -90 to 90 degrees
  std::vector<double> limitsLat = {-90.0, 90.0};
  std::vector<HtlmRegularizationPartParameters> allPartParameters = *params_.parts.value();
  for (auto & partParameters : allPartParameters) {
    HtlmRegularizationPart part(partParameters, allVariables, allLevels, limitsLon, limitsLat);
    applyPart(part, templateFieldSet);
    // Note: this loop applies parts in the order they are listed in the configuration, which will
    // result in overwriting of values if any parts overlap in variable-region space
  }
  for (auto & templateField : templateFieldSet) {
    regularizationFieldSet_.add(templateField);
  }
}

//--------------------------------------------------------------------------------------------------

const double & HtlmRegularizationComponentDependent::getRegularizationValue(
                                                                      const std::string & variable,
                                                                      const size_t & locationN,
                                                                      const size_t & levelN) const {
  auto regularizationFieldArray =
    atlas::array::make_view<double, 2>(regularizationFieldSet_.field(variable));
  return regularizationFieldArray(locationN, levelN);
}

//--------------------------------------------------------------------------------------------------

const bool HtlmRegularizationComponentDependent::AIsInRangeOfB(const double & A,
                                                               const std::vector<double> & B) const
{
  return !((*std::min_element(B.begin(), B.end()) > A)
           || (A > *std::max_element(B.begin(), B.end())));
}

//--------------------------------------------------------------------------------------------------

void HtlmRegularizationComponentDependent::applyPart(const HtlmRegularizationPart & part,
                                                     atlas::FieldSet & templateFieldSet) {
  for (auto & templateField : templateFieldSet) {
    if (!AIsInB(templateField.name(), part.getVariables())) continue;
    auto templateFieldArray = atlas::array::make_view<double, 2>(templateField);
    if (part.containsAllGridPoints()) {
      templateFieldArray.assign(part.getValue());
      continue;
    }
    const auto lonLats = atlas::array::make_view<double, 2>(templateField.functionspace().lonlat());
    for (size_t levelN = 0; levelN < nLevels_; levelN++) {
      if (!AIsInB(levelN, part.getLevels())) continue;
      for (size_t locationN = 0; locationN < nLocations_; locationN++) {
        if (!AIsInRangeOfB(lonLats(locationN, 0), part.getBoundingLons())) continue;
        if (!AIsInRangeOfB(lonLats(locationN, 1), part.getBoundingLats())) continue;
        templateFieldArray(locationN, levelN) = part.getValue();
      }
    }
  }
}

}  // namespace oops
