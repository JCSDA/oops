
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
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// \brief Configuration parameters for a single "part" when constructing regularizationFieldSet_
/// using "parts". Each part specifies which variables and region to apply a particular value to.
class HtlmRegularizationPartParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HtlmRegularizationPartParameters, Parameters);

 public:
  RequiredParameter<double> value{"value", this};
  OptionalParameter<std::vector<std::string>> variables{"variables", this};  // Must be a subset of
  // variables in FieldSet given as argument to HtlmRegularization constructor. Defaults to all.
  OptionalParameter<std::pair<double, double>> boundingLons{"bounding lons", this};  // Band
  // proceeds eastward from lower to higher value. Defaults to including all longitudes.
  OptionalParameter<std::pair<double, double>> boundingLats{"bounding lats", this};  // Defaults to
  // including all latitudes.
  OptionalParameter<std::vector<size_t>> levels{"levels", this};  // Defaults to all.
};

/// \brief A single part for construction by parts.
class HtlmRegularizationPart {
 public:
  HtlmRegularizationPart(const HtlmRegularizationPartParameters &,
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

  const HtlmRegularizationPartParameters params_;
  const double value_;
  const std::vector<std::string> variables_;
  std::vector<size_t> levels_;
  std::pair<double, double> boundingLons_;
  std::pair<double, double> boundingLats_;
  bool containsAllGridPoints_;
};

/// \brief Configuration parameters for the HtlmRegularization class. If just a base value is
/// specified, then this is returned with each call to ::getRegularizationValue(). If a different
/// method, such as by parts (see above), is requested by way of specifying an optional parameter,
/// then an atlas::FieldSet is used for storage of values, and they are retrieved from this.
class HtlmRegularizationParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HtlmRegularizationParameters, Parameters);

 public:
  Parameter<double> baseValue{"base value", 0.0, this};
  OptionalParameter<std::vector<HtlmRegularizationPartParameters>> parts{"parts", this};
};

/// \brief Classes for setting and storage of values used during regularization (ridge regression)
/// when using a hybrid tangent linear model (H-TLM). See https://doi.org/10.1175/MWR-D-20-0088.1
/// Section 4b. Base class handles single-value-only case, derived class handles component-dependent
/// values case.
class HtlmRegularization {
 public:
  explicit HtlmRegularization(const HtlmRegularizationParameters & params)
    : params_(params), baseValue_(params_.baseValue.value()) {}
  virtual ~HtlmRegularization() = default;
  static const std::string classname() {return "oops::HtlmRegularization";}
  virtual const double & getRegularizationValue(const std::string &,
                                                const size_t,
                                                const size_t) const {return baseValue_;}

 protected:
  const HtlmRegularizationParameters params_;
  const double baseValue_;
};

class HtlmRegularizationComponentDependent : public HtlmRegularization {
 public:
  HtlmRegularizationComponentDependent(const HtlmRegularizationParameters &, atlas::FieldSet);
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
