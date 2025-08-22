
// (C) Crown Copyright 2023 Met Office.

// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#ifndef OOPS_GENERIC_HTLMREGULARIZATION_H_
#define OOPS_GENERIC_HTLMREGULARIZATION_H_

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field/FieldSet.h"

namespace oops {

/*
 * Configuration options for a single "part" in regularizationFieldSet_:
 *
 * Each part defines a region to which a regularization value is applied.
 * These are typically specified in a list under a top-level "parts" key.
 *
 * Keys for each part:
 * ─────────────────────────────────────────────────────────────────────────────
 * "value"           : (Required) The regularization value to apply in the specified region.
 *
 * "variables"       : (Optional) List of variable names to apply the value to.
 *                     - Must be a subset of the variables in the FieldSet passed to the
 *                       HtlmRegularization constructor.
 *                     - Defaults to all variables if not specified.
 *
 * "bounding lons"   : (Optional) Longitude bounds as a pair (minLon, maxLon).
 *                     - Defines an eastward band from minLon to maxLon.
 *                     - Defaults to covering all longitudes.
 *
 * "bounding lats"   : (Optional) Latitude bounds as a pair (minLat, maxLat).
 *                     - Defines the north-south extent of the region.
 *                     - Defaults to covering all latitudes.
 *
 * "levels"          : (Optional) List of vertical levels (indices) to apply the value to.
 *                     - Defaults to all levels if not specified.
 */


/// \brief A single part for construction by parts.
class HtlmRegularizationPart {
 public:
  HtlmRegularizationPart(const eckit::Configuration &,
                         const std::vector<std::string> &,
                         const std::vector<size_t> &);
  const std::vector<std::string> & getVariables() const {return variables_;}
  const std::vector<size_t> & getLevels() const {return levels_;}
  const std::pair<double, double> & getBoundingLons() const {return boundingLons_;}
  const std::pair<double, double> & getBoundingLats() const {return boundingLats_;}
  const double & getValue() const {return value_;}
  const bool & containsAllGridPoints() const {return containsAllGridPoints_;}

 private:
  template <typename T>
  bool AIsSubsetOfB(std::vector<T>, std::vector<T>) const;
  bool allOfAAreInRangeOfB(const std::pair<double, double> &,
                           const std::pair<double, double> &) const;

  static const std::pair<double, double> limitsLon;
  static const std::pair<double, double> limitsLat;

  const eckit::LocalConfiguration regPartConfig_;
  const double value_;
  const std::vector<std::string> variables_;
  std::vector<size_t> levels_;
  std::pair<double, double> boundingLons_;
  std::pair<double, double> boundingLats_;
  bool containsAllGridPoints_;
};

/*
 * Configuration options for HtlmRegularization:
 *
 * This configuration controls how regularization values are applied. There are two modes:
 * 1. Uniform regularization using a single "base value".
 * 2. Spatially varying regularization using a list of "parts" (see part-level config above).
 *
 * Keys:
 * ─────────────────────────────────────────────────────────────────────────────
 * "base value" : (Optional) A scalar value used uniformly across all variables and regions.
 *                - Default: 0.0
 *                - If no "parts" are specified, this value is returned by ::getRegularizationValue().
 *
 * "parts"      : (Optional) A list of configuration blocks, each defining a localized regularization.
 */

/// \brief Classes for setting and storage of values used during regularization (ridge regression)
/// when using a hybrid tangent linear model (H-TLM). See https://doi.org/10.1175/MWR-D-20-0088.1
/// Section 4b. Base class handles single-value-only case, derived class handles component-dependent
/// values case.

class HtlmRegularization {
 public:
  explicit HtlmRegularization(const eckit::Configuration & config)
    : regConfig_(config), baseValue_(regConfig_.has("base value") ?
        regConfig_.getDouble("base value") : static_cast<double>(0.0)) {}
  virtual ~HtlmRegularization() = default;
  static const std::string classname() {return "oops::HtlmRegularization";}
  virtual const double & getRegularizationValue(const std::string &,
                                                const size_t,
                                                const size_t) const {return baseValue_;}

 protected:
  const eckit::LocalConfiguration regConfig_;
  const double baseValue_;
};

class HtlmRegularizationComponentDependent : public HtlmRegularization {
 public:
  HtlmRegularizationComponentDependent(const eckit::Configuration &, atlas::FieldSet);
  virtual ~HtlmRegularizationComponentDependent() = default;
  virtual const double & getRegularizationValue(const std::string &,
                                                const size_t,
                                                const size_t) const;

 private:
  template <typename T>
  bool AIsInB(const T &, const std::vector<T> &) const;
  bool AIsInRangeOfB(const double, const std::pair<double, double> &) const;
  void applyPart(const HtlmRegularizationPart &, atlas::FieldSet &);

  const size_t nLevels_;
  const size_t nLocations_;
  atlas::FieldSet regularizationFieldSet_;
};

//--------------------------------------------------------------------------------------------------

template <typename T>
bool HtlmRegularizationPart::AIsSubsetOfB(std::vector<T> A, std::vector<T> B) const {
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  return std::includes(B.begin(), B.end(), A.begin(), A.end());
}

//--------------------------------------------------------------------------------------------------

template <typename T>
bool HtlmRegularizationComponentDependent::AIsInB(const T & A,
                                                        const std::vector<T> & B) const {
  return std::find(B.begin(), B.end(), A) != B.end();
}

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMREGULARIZATION_H_
