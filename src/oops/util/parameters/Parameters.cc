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
#include "oops/util/Timer.h"

#ifdef OOPS_HAVE_NLOHMANN_JSON_SCHEMA_VALIDATOR
#include <boost/algorithm/string/replace.hpp>
#if USE_BOOST_REGEX
#include <boost/regex.hpp>
#define REGEX_NAMESPACE boost
#else
#include <regex>
#define REGEX_NAMESPACE std
#endif
#include "eckit/log/JSON.h"
#include <nlohmann/json-schema.hpp>
#endif

namespace oops {

namespace {

#ifdef OOPS_HAVE_NLOHMANN_JSON_SCHEMA_VALIDATOR

/// \brief JSON schema validation error handler.
///
/// Throws an exception like the default error handler, but tries to make the error message more
/// readable.
class ValidationErrorHandler : public nlohmann::json_schema::basic_error_handler {
 private:
  void error(const nlohmann::json_pointer<nlohmann::basic_json<>> &pointer,
             const nlohmann::json &instance,
             const std::string &message) override {
    // The base class implementation sets an internal flag indicating validation failure
    basic_error_handler::error(pointer, instance, message);

    const std::string editedMessage = postprocessErrorMessage(message);
    const std::string errorLocation = getPath(pointer);

    std::stringstream exceptionMessage;
    exceptionMessage << "Error: YAML validation failed."
                     << "\n  Location:      " << errorLocation
                     << "\n  Invalid value: " << instance.dump()
                     << "\n  Cause:         " << editedMessage;
    throw std::invalid_argument(exceptionMessage.str());
  }

  /// Try to make the error message more readable
  std::string postprocessErrorMessage(const std::string &message) {
    std::string editedMessage = message;
    boost::algorithm::replace_all(editedMessage,
                                  "at least one subschema has failed, but all of them "
                                  "are required to validate - ",
                                  "");

    const REGEX_NAMESPACE::regex additionalPropertyRegex(
                                             "validation failed for additional property '(.*?)': "
                                             "instance invalid as per false-schema");
    editedMessage = REGEX_NAMESPACE::regex_replace(editedMessage, additionalPropertyRegex,
                                       std::string("additional properties are not allowed "
                                       "('$1' was unexpected)"));

    boost::algorithm::replace_all(editedMessage, "instance not found in required enum",
                                  "unrecognized enum value");

    boost::algorithm::replace_all(editedMessage, "instance", "value");

    return editedMessage;
  }

  /// Return the path to the given JSON node
  std::string getPath(const nlohmann::json_pointer<nlohmann::basic_json<>> &pointer) {
    std::string path = pointer.to_string();
    if (path.empty()) {
      // The json_pointer::to_string() method returns an empty string if the pointer refers to
      // the root of the JSON hierarchy. Use a slash instead.
      path = "/";
    }
    return path;
  }
};

void checkStringFormat(const std::string &format, const std::string &value) {
  if (format == "duration") {
    // Matches ISO 8601 representations of durations, for example P1Y (a duration of 1 year),
    // PT1H30M (a duration of 1 h and 30 min), P1YT1S (a duration of 1 year and 1 s)
    static const REGEX_NAMESPACE::regex regex(
          R"(^P(?!$)(\d+Y)?(\d+M)?(\d+W)?(\d+D)?(T(?=\d+[HMS])(\d+H)?(\d+M)?(\d+S)?)?$)");
    REGEX_NAMESPACE::smatch matches;
    if (!REGEX_NAMESPACE::regex_match(value, matches, regex)) {
      throw std::invalid_argument(value + " is not a duration string.");
    }
  } else {
    nlohmann::json_schema::default_string_format_check(format, value);
  }
}
#endif  // OOPS_HAVE_NLOHMANN_JSON_SCHEMA_VALIDATOR

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
  util::Timer timer("oops::Parameters", "deserialize");
  util::CompositePath path;
  deserialize(path, config);
}

void Parameters::deserialize(util::CompositePath &path, const eckit::Configuration &config) {
  for (ParameterBase* child : children_)
    child->deserialize(path, config);
}

void Parameters::serialize(eckit::LocalConfiguration &config) const {
  util::Timer timer("oops::Parameters", "serialize");
  for (const ParameterBase* child : children_) {
    child->serialize(config);
  }
}

eckit::LocalConfiguration Parameters::toConfiguration() const {
  eckit::LocalConfiguration config;
  serialize(config);
  return config;
}

void Parameters::validate(const eckit::Configuration &config) {
#ifdef OOPS_HAVE_NLOHMANN_JSON_SCHEMA_VALIDATOR
  util::Timer timer("oops::Parameters", "validate");
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

  ValidationErrorHandler errorHandler;
  validator.validate(jsonConfig, errorHandler);
#endif
}

void Parameters::validateAndDeserialize(const eckit::Configuration &config) {
  validate(config);
  deserialize(config);
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

  // Allow keys starting with an underscore to be set to arbitrary values. Such keys are useful
  // for defining anchors reused in multiple places in a YAML file.
  schema.extendPatternPropertySchema("^_", PropertyJsonSchema());

  return schema;
}

}  // namespace oops
