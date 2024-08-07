/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <unordered_map>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/utils/Hash.h"

#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
// Note: some of these will not be used and are currently specified for
// testing purposes until the set of staggers is defined.
// If altering the set of staggers, make sure to update the map below and possibly the
// staggerToString function in Variable.cc
enum class VerticalStagger {CENTER,
                            CENTER_WITH_TOP,
                            INTERFACE,
                            TOP_INTERFACE,
                            BOTTOM_INTERFACE};

/// Default vertical stagger (no change to basic name)
static const VerticalStagger defaultVerticalStagger = VerticalStagger::CENTER;
/// Map between all non-default vertical staggers and strings in variable names
static const std::unordered_map<VerticalStagger, std::string> verticalStaggerNameMap {
        {VerticalStagger::INTERFACE,        "_at_interfaces"},
        {VerticalStagger::TOP_INTERFACE,    "_at_top_interfaces"},
        {VerticalStagger::BOTTOM_INTERFACE, "_at_bottom_interfaces"},
        {VerticalStagger::CENTER_WITH_TOP,  "_extended_up_by_1"}};

//-------------------------------------------------------------------------------------
// Enum type for model variable data types
// If altering the set of data types, make sure to update the dataTypeToString function in
// Variable.cc
enum class ModelDataType {
    Byte,
    Int32,
    Int64,
    Real32,
    Real64,
    UInt64
};

static const ModelDataType defaultDataType = ModelDataType::Real64;

//-------------------------------------------------------------------------------------
// Enum type for model variable domain
// If altering the set of domain types, make sure to update the variableDomainToString function in
// Variable.cc
enum class ModelVariableDomain {
    Atmosphere,
    Ocean,
    Land};

static const ModelVariableDomain defaultVariableDomain = ModelVariableDomain::Atmosphere;

// -----------------------------------------------------------------------------
class VariableMetaData : public util::Printable {
 public:
  static std::string staggerToString(const VerticalStagger & stagger);
  static std::string dataTypeToString(const ModelDataType & type);
  static std::string variableDomainToString(const ModelVariableDomain & domain);

  explicit VariableMetaData(const VariableMetaData &) = default;
  explicit VariableMetaData(VariableMetaData &&) = default;
  VariableMetaData(const VerticalStagger & stagger = defaultVerticalStagger,
                   const ModelDataType & type = defaultDataType,
                   const ModelVariableDomain & domain = defaultVariableDomain);
  explicit VariableMetaData(const std::string &);
  VariableMetaData& operator=(const VariableMetaData &) = default;
  VariableMetaData& operator=(VariableMetaData &&) = default;

  bool operator==(const VariableMetaData &) const;
  bool operator!=(const VariableMetaData &) const;

  void print(std::ostream &) const override;

  const VerticalStagger & stagger() const {return stagger_;}
  const ModelDataType & dataType() const {return dataType_;}
  const ModelVariableDomain & domain() const {return domain_;}

 private:
  VerticalStagger stagger_;
  ModelDataType dataType_;
  ModelVariableDomain domain_;
};

// -----------------------------------------------------------------------------
class Variable : public util::Printable {
 public:
  explicit Variable(const std::string &);
  Variable(const std::string &, const VariableMetaData &, const int & levels = -1);
  Variable(const std::string &, const eckit::Configuration &);

  Variable(const Variable &) = default;
  Variable(Variable &&) = default;
  Variable& operator=(const Variable &) = default;
  Variable& operator=(Variable &&) = default;

  // Compares variable names and metadata
  bool operator==(const Variable &) const;
  // Compares variable names and metadata
  bool operator!=(const Variable &) const;
  // Compares only variable names
  bool operator<(const Variable &) const;

  const std::string & name() const {return varName_;}

  int getLevels() const {return levels_;}
  void setLevels(const int & levels) {levels_ = levels;}

  const VariableMetaData & metaData() const {return varMetaData_;}
  const VerticalStagger & stagger() const {return varMetaData_.stagger();}
  const ModelDataType & dataType() const {return varMetaData_.dataType();}

 private:
  void print(std::ostream &) const override;

  std::string varName_;
  VariableMetaData varMetaData_;
  int levels_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

// Specialization of std::hash for Variable
namespace std {
  template <>
  struct hash<oops::Variable> {
    size_t operator()(const oops::Variable & var) const {
      size_t hashName = hash<std::string>{}(var.name());
      return hashName;
    }
  };
}  // namespace std
