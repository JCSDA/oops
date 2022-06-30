/*
 * (C) Copyright 2020 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_REQUIREDPARAMETER_H_
#define OOPS_UTIL_PARAMETERS_REQUIREDPARAMETER_H_

#include <memory>
#include <set>
#include <string>
#include <tuple>
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

/// \brief A mandatory parameter.
///
/// An exception is thrown if the value of this parameter is not found in the Configuration object
/// from which parameters are deserialized.
///
/// \tparam T
///   Type of the value stored in the parameter.
template <typename T>
class RequiredParameter : public ParameterBase {
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
  explicit RequiredParameter(
      const char *name, Parameters *parent,
      std::vector<std::shared_ptr<const ParameterConstraint<T>>> constraints = {})
    : RequiredParameter(name, "", parent, constraints)
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
  explicit RequiredParameter(
      const char *name, const char *description, Parameters *parent,
      std::vector<std::shared_ptr<const ParameterConstraint<T>>> constraints = {})
    : ParameterBase(parent), name_(name), description_(description),
      constraints_(std::move(constraints))
  {}

  /// \brief Constructor with the arguments packaged as a tuple.
  /// This ctor can be called when creating a tuple of RequiredParameters.
  /// \param packagedArgs
  ///   Tuple containing ctor arguments; the tuple elements must exactly match
  ///   (types and ordering) the arguments of the ctor to call.
  template <typename... Args>
  explicit RequiredParameter(const std::tuple<Args...> & packagedArgs)
    : RequiredParameter(packagedArgs, std::index_sequence_for<Args...>{})
  {}

  /// \brief Load the value of this parameter from \p config.
  ///
  /// An exception is thrown if \p config does not contain a key with the name specified in this
  /// parameter's constructor.
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  void serialize(eckit::LocalConfiguration &config) const override;

  ObjectJsonSchema jsonSchema() const override;

  /// \brief The value stored in this parameter.
  ///
  /// An exception is thrown if the value hasn't been loaded from a Configuration yet.
  const T &value() const { return value_.value(); }

  /// \brief The value stored in this parameter.
  ///
  /// An exception is thrown if the value hasn't been loaded from a Configuration yet.
  operator const T &() const { return value_.value(); }

 private:
  std::string name_;
  std::string description_;
  // The value is stored in a boost::optional object because T may not be
  // default-constructible.
  boost::optional<T> value_;
  std::vector<std::shared_ptr<const ParameterConstraint<T>>> constraints_;

  /// \brief Constructor with the arguments packaged as a tuple, taking
  /// index_sequence on input.
  template <typename Tuple, size_t... Is>
  explicit RequiredParameter(const Tuple & packagedArgs, std::index_sequence<Is...> /*meta*/)
      : RequiredParameter(std::get<Is>(packagedArgs)...) {}
};

template <typename T>
void RequiredParameter<T>::deserialize(util::CompositePath &path,
                                       const eckit::Configuration &config) {
  boost::optional<T> newValue = ParameterTraits<T>::get(path, config, name_);
  if (newValue == boost::none)
    throw eckit::BadParameter(path.path() + ": Mandatory parameter '" + name_ + "' not found",
                              Here());
  util::PathComponent component(path, name_);
  for (const std::shared_ptr<const ParameterConstraint<T>> &constraint : constraints_)
    constraint->checkValue(path.path(), newValue.get());  // will throw if constraint is violated
  value_ = std::move(newValue);
}

template <typename T>
void RequiredParameter<T>::serialize(eckit::LocalConfiguration &config) const {
  if (value_ != boost::none)
    ParameterTraits<T>::set(config, name_, *value_);
}

template <typename T>
ObjectJsonSchema RequiredParameter<T>::jsonSchema() const {
  ObjectJsonSchema schema = ParameterTraits<T>::jsonSchema(name_);
  if (description_ != "") {
    schema.extendPropertySchema(name_, {{"description", "\"" + description_ + "\""}});
  }
  schema.require(name_);
  for (const std::shared_ptr<const ParameterConstraint<T>> &constraint : constraints_) {
    PropertyJsonSchema constraintSchema = constraint->jsonSchema();
    schema.extendPropertySchema(name_, std::move(constraintSchema));
  }
  return schema;
}

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_REQUIREDPARAMETER_H_
