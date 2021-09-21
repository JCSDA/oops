/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_OPTIONALPARAMETER_H_
#define OOPS_UTIL_PARAMETERS_OPTIONALPARAMETER_H_

#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/optional.hpp>

#include "eckit/config/Configuration.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ObjectJsonSchema.h"
#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/ParameterConstraint.h"
#include "oops/util/parameters/ParameterTraits.h"

namespace oops {

/// \brief An optional parameter without a default value.
///
/// The interface of this class makes it possible to determine whether the value of this parameter
/// was found in the Configuration object from which parameters were deserialized.
///
/// \tparam T
///   Type of the value stored in the parameter.
template <typename T>
class OptionalParameter : public ParameterBase {
 public:
  /// \brief Constructor.
  ///
  /// \param name
  ///   Name of the key from which this parameter's value will be loaded when parameters are
  ///   deserialized from a Configuration object. Similarly, name of the key to which this
  ///   parameter's value will be saved when parameters are serialized to a Configuration object.
  /// \param parent
  ///   Pointer to the Parameters object representing the collection of options located at
  ///   the same level of the configuration tree as \p name. A call to deserialize() or serialize()
  ///   on that object will automatically trigger a call to deserialize() or serialize() on this
  ///   parameter.
  /// \param constraints
  ///   Zero or more constraints that must be satisfied by the value of this parameter loaded from
  ///   a Configuration object; if that's not the case, an exception will be thrown during
  ///   deserialization.
  explicit OptionalParameter(
      const char *name, Parameters *parent,
      std::vector<std::shared_ptr<const ParameterConstraint<T>>> constraints = {})
    : OptionalParameter(name, "", parent, constraints)
  {}

  /// \brief Constructor.
  ///
  /// \param name
  ///   Name of the key from which this parameter's value will be loaded when parameters are
  ///   deserialized from a Configuration object. Similarly, name of the key to which this
  ///   parameter's value will be saved when parameters are serialized to a Configuration object.
  /// \param description
  ///   Long description of this parameter.
  /// \param parent
  ///   Pointer to the Parameters object representing the collection of options located at
  ///   the same level of the configuration tree as \p name. A call to deserialize() or serialize()
  ///   on that object will automatically trigger a call to deserialize() or serialize() on this
  ///   parameter.
  /// \param constraints
  ///   Zero or more constraints that must be satisfied by the value of this parameter loaded from
  ///   a Configuration object; if that's not the case, an exception will be thrown during
  ///   deserialization.
  explicit OptionalParameter(
      const char *name, const char *description, Parameters *parent,
      std::vector<std::shared_ptr<const ParameterConstraint<T>>> constraints = {})
    : ParameterBase(parent), name_(name), description_(description),
      constraints_(std::move(constraints))
  {}

  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  void serialize(eckit::LocalConfiguration &config) const override;

  ObjectJsonSchema jsonSchema() const override;

  /// \brief The value stored in this parameter, or boost::none if no value is stored.
  const boost::optional<T> &value() const { return value_; }

  /// \brief The value stored in this parameter, or boost::none if no value is stored.
  operator const boost::optional<T> &() const { return value_; }

 private:
  std::string name_;
  std::string description_;
  boost::optional<T> value_;
  std::vector<std::shared_ptr<const ParameterConstraint<T>>> constraints_;
};

template <typename T>
void OptionalParameter<T>::deserialize(util::CompositePath &path,
                                       const eckit::Configuration &config) {
  boost::optional<T> newValue = ParameterTraits<T>::get(path, config, name_);
  if (newValue != boost::none) {
    util::PathComponent component(path, name_);
    for (const std::shared_ptr<const ParameterConstraint<T>> &constraint : constraints_)
      constraint->checkValue(path.path(), newValue.get());  // will throw if constraint is violated
    value_ = std::move(newValue);
  }
}

template <typename T>
void OptionalParameter<T>::serialize(eckit::LocalConfiguration &config) const {
  if (value_ != boost::none)
    ParameterTraits<T>::set(config, name_, *value_);
}

template <typename T>
ObjectJsonSchema OptionalParameter<T>::jsonSchema() const {
  ObjectJsonSchema schema = ParameterTraits<T>::jsonSchema(name_);
  if (description_ != "") {
    schema.extendPropertySchema(name_, {{"description", "\"" + description_ + "\""}});
  }
  for (const std::shared_ptr<const ParameterConstraint<T>> &constraint : constraints_) {
    PropertyJsonSchema constraintSchema = constraint->jsonSchema();
    schema.extendPropertySchema(name_, std::move(constraintSchema));
  }
  return schema;
}


/// \brief Specialisation for `void`.
///
/// OptionalParameter<void> represents YAML options that don't take any values -- only their
/// presence matters. For example, the presence of the following option
///
///     is_defined:
///
/// or equivalently
///
///     is_defined: null
///
/// can be detected by declaring the following member variable of a Parameters subclass:
///
///     oops::OptionalParameter<void> isDefined{"is_defined", this};
///
/// After deserialization, `isDefined.value()` will return true if this option is present in the
/// YAML file and false otherwise.
template <>
class OptionalParameter<void> : public ParameterBase {
 public:
  /// \brief Constructor.
  ///
  /// \param name
  ///   Name of the key from which this parameter's value will be loaded when parameters are
  ///   deserialized from a Configuration object. Similarly, name of the key to which this
  ///   parameter's value will be saved when parameters are serialized to a Configuration object.
  /// \param parent
  ///   Pointer to the Parameters object representing the collection of options located at
  ///   the same level of the configuration tree as \p name. A call to deserialize() or serialize()
  ///   on that object will automatically trigger a call to deserialize() or serialize() on this
  ///   parameter.
  OptionalParameter(const char *name, Parameters *parent)
    : OptionalParameter(name, "", parent)
  {}

  /// \brief Constructor.
  ///
  /// \param name
  ///   Name of the key from which this parameter's value will be loaded when parameters are
  ///   deserialized from a Configuration object. Similarly, name of the key to which this
  ///   parameter's value will be saved when parameters are serialized to a Configuration object.
  /// \param description
  ///   Long description of this parameter.
  /// \param parent
  ///   Pointer to the Parameters object representing the collection of options located at
  ///   the same level of the configuration tree as \p name. A call to deserialize() or serialize()
  ///   on that object will automatically trigger a call to deserialize() or serialize() on this
  ///   parameter.
  OptionalParameter(
      const char *name, const char *description, Parameters *parent)
    : ParameterBase(parent), name_(name), description_(description), value_(false)
  {}

  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  void serialize(eckit::LocalConfiguration &config) const override;

  ObjectJsonSchema jsonSchema() const override;

  /// \brief True if this parameter was present in the deserialized Configuration object, false
  /// otherwise.
  bool value() const { return value_; }

  /// \brief True if this parameter was present in the deserialized Configuration object, false
  /// otherwise.
  operator const bool &() const { return value_; }

 private:
  std::string name_;
  std::string description_;
  bool value_;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_OPTIONALPARAMETER_H_
