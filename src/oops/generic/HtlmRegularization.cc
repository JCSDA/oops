
// (C) Crown Copyright 2023 Met Office.

// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#include "oops/generic/HtlmRegularization.h"

#include <numeric>

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/abor1_cpp.h"

namespace oops {

// Default HtlmRegularizationPart::boundingLons_; 0.0 to 360.0 degrees
const std::pair<double, double> HtlmRegularizationPart::limitsLon = {0.0, 360.0};
// Default HtlmRegularizationPart::boundingLats_; -90.0 to 90.0 degrees
const std::pair<double, double> HtlmRegularizationPart::limitsLat = {-90.0, 90.0};

HtlmRegularizationPart::HtlmRegularizationPart(const HtlmRegularizationPartParameters & params,
                                               const std::vector<std::string> & allVariables,
                                               const std::vector<size_t> & allLevels)
: params_(params), value_(params_.value.value()),
  variables_((params_.variables.value() != boost::none) ?
    *params_.variables.value() : allVariables), levels_(allLevels), boundingLons_(limitsLon),
  boundingLats_(limitsLat), containsAllGridPoints_(true) {
  if (!AIsSubsetOfB(variables_, allVariables)) {
    ABORT("HtlmRegularizationPart: \"variables\" must be a subset of H-TLM update variables");
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
    if (!allOfAAreInRangeOfB(boundingLons_, limitsLon)) {
      ABORT("HtlmRegularizationPart: \"bounding lons\" must be between 0 and 360 degrees");
    }
    containsAllGridPoints_ = false;
  }
  if (params_.boundingLats.value() != boost::none) {
    boundingLats_ = *params.boundingLats.value();
    if (!allOfAAreInRangeOfB(boundingLats_, limitsLat)) {
      ABORT("HtlmRegularizationPart: \"bounding lats\" must be between -90 and 90 degrees");
    }
    containsAllGridPoints_ = false;
  }
}

//--------------------------------------------------------------------------------------------------

bool HtlmRegularizationPart::allOfAAreInRangeOfB(const std::pair<double, double> & A,
                                                 const std::pair<double, double> & B) const {
  return (A.first >= B.first && A.first <= B.second)
          && (A.second >= B.first && A.second <= B.second);
}

//--------------------------------------------------------------------------------------------------

HtlmRegularizationComponentDependent::HtlmRegularizationComponentDependent(
                                                        const HtlmRegularizationParameters & params,
                                                        atlas::FieldSet templateFieldSet)
: HtlmRegularization(params), nLevels_(templateFieldSet[0].shape(1)),
  nLocations_(templateFieldSet[0].shape(0)) {
  for (auto & templateField : templateFieldSet) {
    auto templateArray = atlas::array::make_view<double, 2>(templateField);
    templateArray.assign(baseValue_);
  }
  // Generate default HtlmRegularizationPart::(member variables)...
  // ...for ::variables_; all of the variables in templateFieldSet
  std::vector<std::string> allVariables = templateFieldSet.field_names();
  // ...for ::levels_; all of the model levels indexed from 0
  std::vector<size_t> allLevels(nLevels_);
  std::iota(std::begin(allLevels), std::end(allLevels), 0);
  std::vector<HtlmRegularizationPartParameters> allPartParameters = *params_.parts.value();
  for (auto & partParameters : allPartParameters) {
    HtlmRegularizationPart part(partParameters, allVariables, allLevels);
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
                                                                       const size_t locationN,
                                                                       const size_t levelN) const {
  auto regularizationArray =
    atlas::array::make_view<double, 2>(regularizationFieldSet_.field(variable));
  return regularizationArray(locationN, levelN);
}

//--------------------------------------------------------------------------------------------------

bool HtlmRegularizationComponentDependent::AIsInRangeOfB(
                                                     const double A,
                                                     const std::pair<double, double> & B) const
{
  return !((std::min(B.first, B.second) > A)
           || (A > std::max(B.first, B.second)));
}

//--------------------------------------------------------------------------------------------------

void HtlmRegularizationComponentDependent::applyPart(const HtlmRegularizationPart & part,
                                                     atlas::FieldSet & templateFieldSet) {
  for (auto & templateField : templateFieldSet) {
    if (!AIsInB(templateField.name(), part.getVariables())) continue;
    auto templateArray = atlas::array::make_view<double, 2>(templateField);
    if (part.containsAllGridPoints()) {
      templateArray.assign(part.getValue());
      continue;
    }
    auto lonLats = atlas::array::make_view<double, 2>(templateField.functionspace().lonlat());
    for (size_t levelN = 0; levelN < nLevels_; levelN++) {
      if (!AIsInB(levelN, part.getLevels())) continue;
      for (size_t locationN = 0; locationN < nLocations_; locationN++) {
        if (lonLats(locationN, 0) < 0.0) lonLats(locationN, 0) += 360.0;
        if (!AIsInRangeOfB(lonLats(locationN, 0), part.getBoundingLons())) continue;
        if (!AIsInRangeOfB(lonLats(locationN, 1), part.getBoundingLats())) continue;
        templateArray(locationN, levelN) = part.getValue();
      }
    }
  }
}

}  // namespace oops
