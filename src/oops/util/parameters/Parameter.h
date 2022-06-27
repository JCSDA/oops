/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETER_H_
#define OOPS_UTIL_PARAMETERS_PARAMETER_H_

#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ObjectJsonSchema.h"
#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/ParameterConstraint.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/TypeTraits.h"

namespace oops {

/// \brief A parameter with a default value.
///
/// \note Before declaring a variable of type Parameter<T> where T is something else than a plain
/// old type (int, float etc.), std::map, std::string, std::vector, util::DateTime, util::Duration
/// or a class derived from Parameters, check if there is a dedicated ParameterTraits*.h header
/// containing the specialization of ParameterTraits for type T. If so, include it before the
/// variable declaration. (The "main" ParameterTraits.h header only contains specializations of
/// ParameterTraits for the types listed above).
///
/// \tparam T
///   Type of the value stored in the parameter.
template <typename T>
class Parameter : public ParameterBase {
 public:
  /// \brief Constructor.
  ///
  /// \param name
  ///   Name of the key from which this parameter's value will be loaded when parameters are
  ///   deserialized from a Configuration object. Similarly, name of the key to which this
  ///   parameter's value will be saved when parameters are serialized to a Configuration object.
  /// \param defaultValue
  ///   Default value of the parameter, used if the Configuration object from which parameters are
  ///   deserialized doesn't contain a key with name \p name.
  /// \param parent
  ///   Pointer to the Parameters object representing the collection of options located at
  ///   the same level of the configuration tree as \p name. A call to deserialize() or serialize()
  ///   on that object will automatically trigger a call to deserialize() or serialize() on this
  ///   parameter.
  /// \param constraints
  ///   Zero or more constraints that must be satisfied by the value of this parameter loaded from
  ///   a Configuration object; if that's not the case, an exception will be thrown during
  ///   deserialization.
  Parameter(const char *name, const T& defaultValue, Parameters *parent,
            std::vector<std::shared_ptr<const ParameterConstraint<T>>> constraints = {})
    : Parameter(name, "", defaultValue, parent, constraints)
  {}

  /// \brief Constructor.
  ///
  /// \param name
  ///   Name of the key from which this parameter's value will be loaded when parameters are
  ///   deserialized from a Configuration object. Similarly, name of the key to which this
  ///   parameter's value will be saved when parameters are serialized to a Configuration object.
  /// \param description
  ///   Long description of this parameter.
  /// \param defaultValue
  ///   Default value of the parameter, used if the Configuration object from which parameters are
  ///   deserialized doesn't contain a key with name \p name.
  /// \param parent
  ///   Pointer to the Parameters object representing the collection of options located at
  ///   the same level of the configuration tree as \p name. A call to deserialize() or serialize()
  ///   on that object will automatically trigger a call to deserialize() or serialize() on this
  ///   parameter.
  /// \param constraints
  ///   Zero or more constraints that must be satisfied by the value of this parameter loaded from
  ///   a Configuration object; if that's not the case, an exception will be thrown during
  ///   deserialization.
  Parameter(const char *name, const char *description, const T& defaultValue, Parameters *parent,
            std::vector<std::shared_ptr<const ParameterConstraint<T>>> constraints = {})
    : ParameterBase(parent), name_(name), description_(description), value_(defaultValue),
      defaultValue_(defaultValue), constraints_(std::move(constraints))
  {}

  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  void serialize(eckit::LocalConfiguration &config) const override;

  ObjectJsonSchema jsonSchema() const override;

  /// \brief The value stored in this parameter.
  const T &value() const { return value_; }

  /// \brief The value stored in this parameter.
  operator const T &() const { return value_; }

 private:
  std::string name_;
  std::string description_;
  T value_;
  T defaultValue_;
  std::vector<std::shared_ptr<const ParameterConstraint<T>>> constraints_;
};

template <typename T>
void Parameter<T>::deserialize(util::CompositePath &path,
                               const eckit::Configuration &config) {
  boost::optional<T> newValue = ParameterTraits<T>::get(path, config, name_);
  if (newValue != boost::none) {
    util::PathComponent component(path, name_);
    for (const std::shared_ptr<const ParameterConstraint<T>> &constraint : constraints_)
      constraint->checkValue(path.path(), newValue.get());  // will throw if constraint is violated
    value_ = std::move(newValue.get());
  }
}

template <typename T>
void Parameter<T>::serialize(eckit::LocalConfiguration &config) const {
  ParameterTraits<T>::set(config, name_, value_);
}

template<class, class = void>
struct HasValueAsJson : std::false_type {};

template<class T>
struct HasValueAsJson<T, cpp17::void_t<
  decltype(ParameterTraits<T>::valueAsJson(std::declval<T>()))>> : std::true_type
{};

// Overload selected when ParameterTraits<T> has a static member function `valueAsJson()`
// taking a parameter of type `T`.
template <typename T>
typename std::enable_if<HasValueAsJson<T>::value, void>::type
maybeExtendPropertySchemaWithDefault(ObjectJsonSchema &schema,
                                     const std::string &propertyName,
                                     const T &defaultValue)
{
  std::string defaultValueAsJson = ParameterTraits<T>::valueAsJson(defaultValue);
  if (defaultValueAsJson != "") {
    schema.extendPropertySchema(propertyName, {{"default", defaultValueAsJson}});
  }
}

// No-op overload selected otherwise
template <typename T>
typename std::enable_if<!HasValueAsJson<T>::value, void>::type
maybeExtendPropertySchemaWithDefault(ObjectJsonSchema &schema,
                                     const std::string &propertyName,
                                     const T &defaultValue)
{}

template <typename T>
ObjectJsonSchema Parameter<T>::jsonSchema() const {
  ObjectJsonSchema schema = ParameterTraits<T>::jsonSchema(name_);
  if (description_ != "") {
    schema.extendPropertySchema(name_, {{"description", "\"" + description_ + "\""}});
  }
  maybeExtendPropertySchemaWithDefault(schema, name_, defaultValue_);
  for (const std::shared_ptr<const ParameterConstraint<T>> &constraint : constraints_) {
    PropertyJsonSchema constraintSchema = constraint->jsonSchema();
    schema.extendPropertySchema(name_, std::move(constraintSchema));
  }
  return schema;
}

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETER_H_
