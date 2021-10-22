/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARAMETERS_OBJECTJSONSCHEMA_H_
#define OOPS_UTIL_PARAMETERS_OBJECTJSONSCHEMA_H_

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "oops/util/parameters/PropertyJsonSchema.h"

namespace oops {

struct ConditionalObjectJsonSchema;

/// \brief A representation of a JSON schema used to validate a JSON node of type "object" (a map
/// in C++ parlance).
///
/// Some keywords defined by the JSON Schema specification, but not currently required to describe
/// the structure of JEDI configuration files (e.g. oneOf, anyOf) are omitted.
class ObjectJsonSchema {
 public:
  typedef std::map<std::string, PropertyJsonSchema> PropertyJsonSchemas;

  /// \brief Create a JSON schema describing a node expected to contain certain properties ("keys").
  ///
  /// \param properties
  ///   Map of property names to JSON schemas used to validate these properties.
  /// \param required
  ///   Names of properties that must be present in a valid JSON node (all others are optional).
  /// \param additionalProperties
  ///   True if a valid node may have properties other than those listed in \p properties,
  ///   false otherwise.
  explicit ObjectJsonSchema(PropertyJsonSchemas properties = {},
                            std::set<std::string> required = {},
                            bool additionalProperties = false);

  /// \brief Create a JSON schema describing a node expected to contain a "special" property whose
  /// value determines the list of other properties expected to be present.
  ///
  /// \param selector
  ///   Name of the property whose value determines which other properties may be present.
  /// \param variants
  ///   `variants[key]` is the JSON schema against which the node will be validated if the selector
  ///   property has the value `value`.
  ObjectJsonSchema(const std::string &selector,
                   std::map<std::string, ObjectJsonSchema> variants);

  /// \brief Create a JSON schema representing a logical conjunction of a number of conditional
  /// JSON schemas.
  explicit ObjectJsonSchema(std::vector<ConditionalObjectJsonSchema> allOf);

  /// \brief Map of property names to JSON schemas used to validate these properties.
  const PropertyJsonSchemas &properties() const { return properties_; }

  /// \brief Map of regular expressions to JSON schemas used to validate properties with names
  /// matching those expressions.
  const PropertyJsonSchemas &patternProperties() const { return patternProperties_; }

  /// \brief Names of properties that must be present in a valid JSON node.
  const std::set<std::string> &required() const { return required_; }

  /// \brief True if a valid JSON node may have properties other than those listed in
  /// \c properties, false otherwise.
  bool additionalProperties() const { return additionalProperties_; }

  /// \brief Extra JSON schemas to which a valid JSON node must conform.
  const std::vector<ConditionalObjectJsonSchema> &allOf() const { return allOf_; }

  /// \brief Return true if this is an empty (trivial) schema.
  bool empty() const;

  /// \brief Return a string containing the JSON schema represented by this object.
  ///
  /// \param includeSchemaKeyword
  ///   If true, the returned schema will contain a \c $schema keyword at the top level. This is
  ///   typically useful only for the top-level schema, and not for subschemas.
  ///   (This feature is currently disabled.)
  std::string toString(bool includeSchemaKeyword = false) const;

  /// \brief Convert this object to a PropertyJsonSchema representing the same schema.
  PropertyJsonSchema toPropertyJsonSchema() const;

  /// \brief Combine this schema with another schema (in place).
  ///
  /// This is done as follows:
  ///
  /// - property schemas defined by the other schema are added to the property schemas defined by
  ///   this schema; if both this schema and the other schema define a schema for a property and
  ///   one of the schemas includes the other, the "larger" (more constraining) schema is selected;
  ///   otherwise the two schemas are deemed to be in conflict and an exception is thrown
  ///
  /// - properties required by the other schema are added to those required by this schema
  ///
  /// - additional properties are allowed if either this or the other schema allows them
  ///
  /// - the list of extra (conditional) JSON schemas to which a valid JSON node must conform is
  ///   effectively the Cartesian product of the lists taken from this schema and the other schema.
  ///   See the implementation for more details.
  void combineWith(const ObjectJsonSchema& other);

  /// \brief Add the key-value pairs from \p schema into the schema used to validate the
  /// property \p property.
  ///
  /// Pre-existing key-value pairs take precedence, i.e. a key-value pair in \p schema will be
  /// ignored if the property schema already contains that key.
  void extendPropertySchema(const std::string &property, const PropertyJsonSchema &schema);

  /// \brief Add the key-value pairs from \p schema into the schema used to validate properties
  /// matching the regular expression \p property.
  ///
  /// Pre-existing key-value pairs take precedence, i.e. a key-value pair in \p schema will be
  /// ignored if the property schema already contains that key.
  void extendPatternPropertySchema(const std::string &property, const PropertyJsonSchema &schema);

  /// \brief Mark property \p property as required.
  void require(const std::string &property);

 private:
  std::string propertiesToString() const;
  std::string patternPropertiesToString() const;
  std::string requiredToString() const;
  std::string additionalPropertiesToString() const;
  std::string allOfToString() const;

  static std::string propertyJsonSchemasToString(const PropertyJsonSchemas &properties);

  void combinePropertiesWith(const ObjectJsonSchema& other);
  void combinePatternPropertiesWith(const ObjectJsonSchema& other);
  void combineRequiredWith(const ObjectJsonSchema& other);
  void combineAdditionalPropertiesWith(const ObjectJsonSchema& other);
  void combineAllOfWith(const ObjectJsonSchema& other);

  static void combinePropertyJsonSchemasWith(PropertyJsonSchemas &properties,
                                             const PropertyJsonSchemas& otherProperties);

 private:
  /// \brief Maps property names to JSON schemas used to validate these properties.
  PropertyJsonSchemas properties_;
  /// \brief Maps regular expressions to JSON schemas used to validate properties whose names match
  /// those regular expressions.
  PropertyJsonSchemas patternProperties_;
  /// \brief Names of properties that must be present in a valid JSON node.
  std::set<std::string> required_;
  /// \brief True if a valid JSON node may have properties other than those listed in
  /// \c properties, false otherwise.
  bool additionalProperties_ = false;
  /// \brief Extra JSON schemas to which a valid JSON node must conform.
  std::vector<ConditionalObjectJsonSchema> allOf_;
};

/// \brief A conditional JSON schema.
///
/// A JSON node conforms to this schema if either
/// - it conforms to the schema defined by the \c if_ member variable and to that defined by the
///   \c then member variable, or
/// - it does not conform to the schema defined by the \c if_ member, but it conforms to the schema
///   defined by the \c else_ member.
struct ConditionalObjectJsonSchema {
  ConditionalObjectJsonSchema() = default;

  /// Note: omitting \p elseSchema leads to a schema with a trivial \c else schema (matching
  /// everything).
  ConditionalObjectJsonSchema(ObjectJsonSchema ifSchema, ObjectJsonSchema thenSchema,
                              ObjectJsonSchema elseSchema = ObjectJsonSchema())
    : if_(std::move(ifSchema)), then(std::move(thenSchema)), else_(std::move(elseSchema))
  {}

  ObjectJsonSchema if_;
  ObjectJsonSchema then;
  ObjectJsonSchema else_;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_OBJECTJSONSCHEMA_H_
