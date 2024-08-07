/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/Variable.h"

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "eckit/types/Types.h"
#include "eckit/utils/Hash.h"

#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace oops {

// -----------------------------------------------------------------------------
VariableMetaData::VariableMetaData(const VerticalStagger & stagger, const ModelDataType & type,
                                   const ModelVariableDomain & domain)
  : stagger_(stagger), dataType_(type), domain_(domain) {}

// -----------------------------------------------------------------------------
VariableMetaData::VariableMetaData(const std::string & name) {
  // Parse the name to get the stagger and data type
  // For now, just set the stagger to default, data type to double, and domain to atmosphere
  stagger_ = defaultVerticalStagger;
  dataType_ = defaultDataType;
  domain_ = defaultVariableDomain;
}

// -----------------------------------------------------------------------------

std::string VariableMetaData::staggerToString(const VerticalStagger & stagger) {
  if (verticalStaggerNameMap.find(stagger) != verticalStaggerNameMap.end()) {
    return verticalStaggerNameMap.at(stagger);
  }
  if (stagger == VerticalStagger::CENTER) return "center";
  return "unknown stagger: " + std::to_string(static_cast<int>(stagger));
}

// -----------------------------------------------------------------------------

std::string VariableMetaData::dataTypeToString(const ModelDataType & type) {
  switch (type) {
    case ModelDataType::Byte: return "byte";
    case ModelDataType::Int32: return "int32";
    case ModelDataType::Int64: return "int64";
    case ModelDataType::Real32: return "real32";
    case ModelDataType::Real64: return "real64";
    case ModelDataType::UInt64: return "uint64";
    return "unknown data type: " + std::to_string(static_cast<int>(type));
  }
}

// -----------------------------------------------------------------------------

std::string VariableMetaData::variableDomainToString(const ModelVariableDomain & domain) {
  switch (domain) {
    case ModelVariableDomain::Atmosphere: return "atmosphere";
    case ModelVariableDomain::Ocean: return "ocean";
    case ModelVariableDomain::Land: return "land";
    return "unknown model variable domain: " + std::to_string(static_cast<int>(domain));
  }
}

// -----------------------------------------------------------------------------

bool VariableMetaData::operator==(const VariableMetaData & rhs) const {
  return (stagger_  == rhs.stagger_ &&
          dataType_ == rhs.dataType_ &&
          domain_   == rhs.domain_);
}

// -----------------------------------------------------------------------------

bool VariableMetaData::operator!=(const VariableMetaData & rhs) const {
  return (!(*this == rhs));
}

// -----------------------------------------------------------------------------

void VariableMetaData::print(std::ostream & os) const {
  os << "Stagger: " << staggerToString(stagger_) << ", "
     << "Data type: " << dataTypeToString(dataType_) << ", "
     << "Domain: " << variableDomainToString(domain_);
}

// -----------------------------------------------------------------------------

Variable::Variable(const std::string & name, const VariableMetaData & md, const int & levels)
  : varName_(name), varMetaData_(md), levels_(levels) {}

// -----------------------------------------------------------------------------
Variable::Variable(const std::string & name)
  : varName_(name), varMetaData_(name), levels_(-1) {}

// -----------------------------------------------------------------------------
Variable::Variable(const std::string & name, const eckit::Configuration & conf)
  : varName_(name), varMetaData_(name) {
    if (conf.has("levels")) {
      levels_ = conf.getInt("levels");
    } else {
      levels_ = -1;
    }
}

// -----------------------------------------------------------------------------

bool Variable::operator==(const Variable & rhs) const {
  // For now, ignore levels in comparison if either side has not been assigned a valid value
  bool levelsEqual = true;
  if (levels_ >= 0 && rhs.levels_ >= 0) {
    levelsEqual = (levels_ == rhs.levels_);
  }
  return (varName_ == rhs.varName_ && varMetaData_ == rhs.varMetaData_ && levelsEqual);
}

// -----------------------------------------------------------------------------

bool Variable::operator!=(const Variable & rhs) const {
  return (!(*this == rhs));
}

// -----------------------------------------------------------------------------

bool Variable::operator<(const Variable & other) const {
  return (varName_ < other.varName_);
}

// -----------------------------------------------------------------------------

void Variable::print(std::ostream & os) const {
  os << varName_;
  // When debugging, maybe additional info desired, and this line can be uncommented
  // os << " (" << varMetaData_ << "), levels: " << levels_;
}

}  // namespace oops
