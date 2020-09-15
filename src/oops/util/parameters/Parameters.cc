/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <iostream>
#include <iterator>
#include <utility>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/value/Value.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ObjectJsonSchema.h"
#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/Parameters.h"

#ifdef OOPS_HAVE_NLOHMANN_JSON_SCHEMA_VALIDATOR
#include <regex>
#include "eckit/log/JSON.h"
#include <nlohmann/json-schema.hpp>
#endif

namespace oops {

namespace {

#ifdef OOPS_HAVE_NLOHMANN_JSON_SCHEMA_VALIDATOR
void checkStringFormat(const std::string &format, const std::string &value) {
  if (format == "duration") {
    // Matches ISO 8601 representations of durations, for example P1Y (a duration of 1 year),
    // PT1H30M (a duration of 1 h and 30 min), P1YT1S (a duration of 1 year and 1 s)
    static const std::regex regex(
          R"(^P(?!$)(\d+Y)?(\d+M)?(\d+W)?(\d+D)?(T(?=\d+[HMS])(\d+H)?(\d+M)?(\d+S)?)?$)");

    std::smatch matches;
    if (!std::regex_match(value, matches, regex)) {
      throw std::invalid_argument(value + " is not a duration string.");
    }
  } else {
    nlohmann::json_schema::default_string_format_check(format, value);
  }
}
#endif

}  // namespace

void Parameters::print(std::ostream &os) const {
  eckit::LocalConfiguration config;
  serialize(config);
  os << config;
}

void Parameters::registerChild(ParameterBase &parameter) {
  children_.push_back(&parameter);
}

void Parameters::deserialize(const eckit::Configuration &config) {
  util::CompositePath path;
  deserialize(path, config);
}

void Parameters::deserialize(util::CompositePath &path, const eckit::Configuration &config) {
  for (ParameterBase* child : children_)
    child->deserialize(path, config);
}

void Parameters::serialize(eckit::LocalConfiguration &config) const {
  for (const ParameterBase* child : children_) {
    child->serialize(config);
  }
}

void Parameters::validate(const eckit::Configuration &config) {
#ifdef OOPS_HAVE_NLOHMANN_JSON_SCHEMA_VALIDATOR
  std::string strSchema = jsonSchema().toString(true /*includeSchemaKeyword?*/);
  nlohmann::json jsonSchema = nlohmann::json::parse(strSchema);

  std::stringstream jsonStream;
  {
    eckit::JSON json(jsonStream);
    json << config;
  }
  std::string jsonString = jsonStream.str();
  if (jsonString == "null") {
    // Stop the validator from complaining the top-level JSON node is not of the correct type
    jsonString = "{}";
  }

  const nlohmann::json jsonConfig = nlohmann::json::parse(jsonString);

  nlohmann::json_schema::json_validator validator(nullptr, checkStringFormat);
  validator.set_root_schema(jsonSchema);
  validator.validate(jsonConfig);
#endif
}

bool Parameters::isValidationSupported() {
#ifdef OOPS_HAVE_NLOHMANN_JSON_SCHEMA_VALIDATOR
  return true;
#else
  return false;
#endif
}

ObjectJsonSchema Parameters::jsonSchema() const {
  std::vector<ObjectJsonSchema> childSchemas;
  childSchemas.reserve(children_.size());

  ObjectJsonSchema schema;
  for (const ParameterBase* child : children_) {
    ObjectJsonSchema childSchema = child->jsonSchema();
    schema.combineWith(childSchema);
  }
  return schema;
}

}  // namespace oops
