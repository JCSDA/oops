
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

HtlmRegularizationPart::HtlmRegularizationPart(const eckit::Configuration & config,
                                               const std::vector<std::string> & allVariables,
                                               const std::vector<size_t> & allLevels)
: regPartConfig_(config), value_(regPartConfig_.getDouble("value")),
  variables_((regPartConfig_.has("variables")) ?
    regPartConfig_.getStringVector("variables") : allVariables),
                       levels_(allLevels), boundingLons_(limitsLon),
  boundingLats_(limitsLat), containsAllGridPoints_(true) {
  if (!AIsSubsetOfB(variables_, allVariables)) {
    ABORT("HtlmRegularizationPart: \"variables\" must be a subset of H-TLM update variables");
  }
  if (regPartConfig_.has("levels")) {
    levels_ = regPartConfig_.getUnsignedVector("levels");
    if (!AIsSubsetOfB(levels_, allLevels)) {
      ABORT("HtlmRegularizationPart: \"levels\" must be a subset of model levels indexed from 0");
    }
    containsAllGridPoints_ = false;
  }
  if (regPartConfig_.has("bounding lons")) {
    std::vector<double> blons = regPartConfig_.getDoubleVector("bounding lons");
    boundingLons_ = {blons[0], blons[1]};
    if (!allOfAAreInRangeOfB(boundingLons_, limitsLon)) {
      ABORT("HtlmRegularizationPart: \"bounding lons\" must be between 0 and 360 degrees");
    }
    containsAllGridPoints_ = false;
  }
  if (regPartConfig_.has("bounding lats")) {
    std::vector<double> blats = regPartConfig_.getDoubleVector("bounding lats");
    boundingLats_ = {blats[0], blats[1]};
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
                                                        const eckit::Configuration & config,
                                                        atlas::FieldSet templateFieldSet)
: HtlmRegularization(config), nLevels_(templateFieldSet[0].shape(1)),
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
  std::vector<eckit::LocalConfiguration> allPartConfigs;
  config.get("parts", allPartConfigs);
  for (auto & partConfig : allPartConfigs) {
    HtlmRegularizationPart part(partConfig, allVariables, allLevels);
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
