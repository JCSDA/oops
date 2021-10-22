/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/parameters/ObjectJsonSchema.h"

#include <algorithm>
#include <sstream>
#include <utility>

#include "eckit/exception/Exceptions.h"
#include "eckit/log/Channel.h"
#include "oops/util/stringFunctions.h"

namespace oops {

namespace {

/// \brief Given two *non-empty* vectors of conditional schemas, returns their Cartesian product.
std::vector<ConditionalObjectJsonSchema> cartesianProduct(
    const std::vector<ConditionalObjectJsonSchema> &schemasA,
    const std::vector<ConditionalObjectJsonSchema> &schemasB)
{
  std::vector<ConditionalObjectJsonSchema> result;
  result.reserve(schemasA.size() * schemasB.size());
  for (const ConditionalObjectJsonSchema &schemaA : schemasA)
    for (const ConditionalObjectJsonSchema &schemaB : schemasB) {
      ConditionalObjectJsonSchema combinedSchema = schemaA;
      combinedSchema.if_.combineWith(schemaB.if_);
      combinedSchema.then.combineWith(schemaB.then);
      combinedSchema.else_.combineWith(schemaB.else_);
      result.push_back(std::move(combinedSchema));
    }
  return result;
}

/// \brief For each key `key` in `source` that doesn't already belong to `target`, set
/// `target[key] = source[key]`.
void combinePropertySchemaWith(PropertyJsonSchema &target, const PropertyJsonSchema& source) {
  for (const std::pair<const std::string, std::string> &keyAndValue : source) {
    auto it = target.find(keyAndValue.first);
    // Proritise existing key-value pairs: insert only if it's a new key
    if (it == target.end()) {
      target.insert(keyAndValue);
    }
  }
}

/// \brief Return `s` enclosed in quotes.
std::string quote(const std::string &s) {
  std::string result = "\"";
  result += s;
  result += "\"";
  return result;
}

/// \brief Return the first element of the pair `keyAndValue` enclosed in quotes.
std::string quotedKey(const std::pair<const std::string, ObjectJsonSchema> &keyAndValue) {
  return quote(keyAndValue.first);
}

}  // namespace

ObjectJsonSchema::ObjectJsonSchema(std::map<std::string, PropertyJsonSchema> properties,
                                   std::set<std::string> required,
                                   bool additionalProperties)
  : properties_(properties), required_(required), additionalProperties_(additionalProperties)
{}

ObjectJsonSchema::ObjectJsonSchema(const std::string &selector,
                                   std::map<std::string, ObjectJsonSchema> variants) {
  PropertyJsonSchema &selectorProperties = properties_[selector];
  selectorProperties["enum"] = '[' +
      util::stringfunctions::join(", ", variants.begin(), variants.end(), quotedKey) +
      ']';
  selectorProperties["type"] = R"(["string", "number"])";

  allOf_.reserve(variants.size());
  for (std::pair<const std::string, ObjectJsonSchema> &keyAndSchema : variants) {
    ObjectJsonSchema if_({{selector, {{"const", quote(keyAndSchema.first)}}}}, {}, true);
    ObjectJsonSchema &then = keyAndSchema.second;
    allOf_.push_back(ConditionalObjectJsonSchema(std::move(if_), std::move(then)));
  }
}

ObjectJsonSchema::ObjectJsonSchema(std::vector<ConditionalObjectJsonSchema> allOf)
  : allOf_(allOf)
{}

bool ObjectJsonSchema::empty() const {
  return properties_.empty() && required_.empty() && allOf_.empty();
}

std::string ObjectJsonSchema::toString(bool /*includeSchemaKeyword*/) const {
  PropertyJsonSchema schema = toPropertyJsonSchema();
  // TODO(wsmigaj): handle the includeSchemaKeyword parameter.
  return oops::toString(schema);
}

PropertyJsonSchema ObjectJsonSchema::toPropertyJsonSchema() const {
  PropertyJsonSchema result;

  result["type"] = R"("object")";

  std::string str = propertiesToString();
  if (!str.empty())
    result["properties"] = str;

  str = patternPropertiesToString();
  if (!str.empty())
    result["patternProperties"] = str;

  str = requiredToString();
  if (!str.empty())
    result["required"] = str;

  str = additionalPropertiesToString();
  if (!str.empty())
    result["additionalProperties"] = str;

  str = allOfToString();
  if (!str.empty())
    result["allOf"] = str;

  return result;
}

std::string ObjectJsonSchema::propertiesToString() const {
  return propertyJsonSchemasToString(properties_);
}

std::string ObjectJsonSchema::patternPropertiesToString() const {
  return propertyJsonSchemasToString(patternProperties_);
}

std::string ObjectJsonSchema::propertyJsonSchemasToString(const PropertyJsonSchemas &properties) {
  std::stringstream stream;
  if (!properties.empty()) {
    eckit::Channel channel;
    channel.setStream(stream);

    channel << "{\n";
    bool needsCommaAndNewline = false;
    {
      eckit::AutoIndent indent(channel);
      for (const auto &nameAndSchema : properties) {
        if (needsCommaAndNewline)
          channel << ",\n";
        channel << '"' << nameAndSchema.first << '"'
                << ": " << oops::toString(nameAndSchema.second);
        needsCommaAndNewline = true;
      }
    }
    if (needsCommaAndNewline)
      channel << '\n';  // this was the last item, so comma is not needed, just a newline
    channel << '}';
  }
  return stream.str();
}

std::string ObjectJsonSchema::requiredToString() const {
  std::stringstream stream;
  if (!required_.empty()) {
    eckit::Channel channel;
    channel.setStream(stream);

    channel << '[';
    bool needsCommaAndSpace = false;
    for (const std::string &name : required_) {
      if (needsCommaAndSpace)
        channel << ", ";
      channel << '"' << name << '"';
      needsCommaAndSpace = true;
    }
    channel << ']';
  }
  return stream.str();
}

std::string ObjectJsonSchema::additionalPropertiesToString() const {
  std::string result;
  // If there is an allOf clause, then additionalProperties are enabled or disabled in the
  // conditional expressions embedded in that clause rather than at the top level.
  if (!additionalProperties_ && allOf_.empty()) {
    result = "false";
  }
  return result;
}

std::string ObjectJsonSchema::allOfToString() const {
  std::stringstream stream;
  if (!allOf_.empty()) {
    eckit::Channel channel;
    channel.setStream(stream);

    channel << "[\n";
    {
      eckit::AutoIndent indent(channel);
      bool needsCommaAndNewline = false;
      for (const auto &conditionalSchema : allOf_) {
        if (needsCommaAndNewline)
          channel << ",\n";
        channel << "{\n";
        {
          eckit::AutoIndent indent(channel);
          channel << R"("if": )";
          channel << conditionalSchema.if_.toString();
          if (!conditionalSchema.then.empty()) {
            channel << ",\n";
            channel << R"("then": )";
            channel << conditionalSchema.then.toString();
          }
          if (!conditionalSchema.else_.empty()) {
            channel << ",\n";
            channel << R"("else": )";
            channel << conditionalSchema.else_.toString();
          }
        }
        channel << "\n}";
        needsCommaAndNewline = true;
      }
    }
    channel << ']';
  }
  return stream.str();
}

void ObjectJsonSchema::combineWith(const ObjectJsonSchema& other) {
  combinePropertiesWith(other);
  combinePatternPropertiesWith(other);
  combineRequiredWith(other);
  combineAdditionalPropertiesWith(other);
  // This call needs to be at the end.
  combineAllOfWith(other);
}

void ObjectJsonSchema::combinePropertiesWith(const ObjectJsonSchema& other) {
  combinePropertyJsonSchemasWith(properties_, other.properties());
}

void ObjectJsonSchema::combinePatternPropertiesWith(const ObjectJsonSchema& other) {
  combinePropertyJsonSchemasWith(patternProperties_, other.patternProperties());
}

void ObjectJsonSchema::combinePropertyJsonSchemasWith(PropertyJsonSchemas &properties,
                                                      const PropertyJsonSchemas& otherProperties) {
  for (const std::pair<const std::string, PropertyJsonSchema> &keyAndValue : otherProperties) {
    const std::string &key = keyAndValue.first;
    auto it = properties.find(key);
    if (it == properties.end()) {
      // This schema doesn't define a schema for the property with this key.
      // Import the definition from the other schema.
      properties.insert(keyAndValue);
    } else {
      // This schema already defines a schema for the property with this key.
      PropertyJsonSchema &ourSchema = it->second;
      const PropertyJsonSchema &otherSchema = keyAndValue.second;

      if (std::includes(ourSchema.begin(), ourSchema.end(),
                        otherSchema.begin(), otherSchema.end())) {
        // Our schema already includes the other schema. Nothing needs to be done.
      } else if (std::includes(otherSchema.begin(), otherSchema.end(),
                               ourSchema.begin(), ourSchema.end())) {
        // The other schema includes our schema. Replace our schema with the other schema.
        ourSchema = otherSchema;
      } else {
        // None of the schemas is a subset of the other.
        ASSERT_MSG(false, "Schemas to be combined define conflicting schemas "
                          "for the property '" + key + "'");
      }
    }
  }
}

void ObjectJsonSchema::combineRequiredWith(const ObjectJsonSchema& other) {
  required_.insert(other.required_.begin(), other.required_.end());
}

void ObjectJsonSchema::combineAdditionalPropertiesWith(const ObjectJsonSchema &other) {
  additionalProperties_ = additionalProperties_ || other.additionalProperties_;
}

void ObjectJsonSchema::combineAllOfWith(const ObjectJsonSchema& other) {
  if (!other.allOf_.empty()) {
    if (allOf_.empty())
      allOf_ = other.allOf_;
    else
      allOf_ = cartesianProduct(allOf_, other.allOf_);
  }

  // Note: this object's properties, required and additionalProperties member variables have
  // already been combined with those of the other schema.
  for (ConditionalObjectJsonSchema &conditionalSchema : allOf_) {
    if (!conditionalSchema.then.empty()) {
      conditionalSchema.then.combinePropertiesWith(*this);
      conditionalSchema.then.combinePatternPropertiesWith(*this);
      conditionalSchema.then.combineRequiredWith(*this);
      conditionalSchema.then.combineAdditionalPropertiesWith(*this);
    }
    if (!conditionalSchema.else_.empty()) {
      conditionalSchema.else_.combinePropertiesWith(*this);
      conditionalSchema.else_.combinePatternPropertiesWith(*this);
      conditionalSchema.else_.combineRequiredWith(*this);
      conditionalSchema.else_.combineAdditionalPropertiesWith(*this);
    }
  }
}

void ObjectJsonSchema::extendPropertySchema(const std::string &property,
                                            const PropertyJsonSchema &schema) {
  combinePropertySchemaWith(properties_[property], schema);
  for (ConditionalObjectJsonSchema &conditionalSchema : allOf_)
    combinePropertySchemaWith(conditionalSchema.then.properties_[property], schema);
}

void ObjectJsonSchema::extendPatternPropertySchema(const std::string &property,
                                                   const PropertyJsonSchema &schema) {
  combinePropertySchemaWith(patternProperties_[property], schema);
  for (ConditionalObjectJsonSchema &conditionalSchema : allOf_)
    combinePropertySchemaWith(conditionalSchema.then.patternProperties_[property], schema);
}


void ObjectJsonSchema::require(const std::string &property) {
  required_.insert(property);
  for (ConditionalObjectJsonSchema &conditionalSchema : allOf_) {
    if (!conditionalSchema.then.empty())
      conditionalSchema.then.required_.insert(property);
    if (!conditionalSchema.else_.empty())
      conditionalSchema.else_.required_.insert(property);
  }
}

}  // namespace oops
