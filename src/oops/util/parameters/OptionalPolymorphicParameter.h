/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_OPTIONALPOLYMORPHICPARAMETER_H_
#define OOPS_UTIL_PARAMETERS_OPTIONALPOLYMORPHICPARAMETER_H_

#include <memory>
#include <set>
#include <string>
#include <utility>

#include <boost/optional.hpp>

#include "oops/util/parameters/ObjectJsonSchema.h"
#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/PolymorphicParameterTraits.h"

namespace util {
class CompositePath;
}

namespace oops {

/// \brief An optional parameter storing an instance of a subclass of `PARAMETERS` (itself a
/// subclass of Parameters).
///
/// The subclass is determined at runtime by loading the value of a specified key from a
/// Configuration object. If that key is not present, the parameter remains empty.
///
/// See RequiredPolymorphicParameter for general information about polymorphic parameters and an
/// explanation of the `FACTORY` template parameter.
///
/// \see RequiredPolymorphicParameter, PolymorphicParameter
template <typename PARAMETERS, typename FACTORY>
class OptionalPolymorphicParameter : public ParameterBase {
 public:
  /// \brief Constructor.
  ///
  /// \param name
  ///   Name of the configuration key whose value determines the concrete subclass of `PARAMETERS`
  ///   created during deserialization.
  /// \param parent
  ///   Pointer to the Parameters object representing the collection of options located at
  ///   the same level of the configuration tree as `name`. A call to deserialize() or serialize()
  ///   on that object will automatically trigger a call to deserialize() or serialize() on this
  ///   parameter.
  explicit OptionalPolymorphicParameter(const char *name, Parameters *parent)
    : OptionalPolymorphicParameter(name, "", parent)
  {}

  /// \brief Constructor.
  ///
  /// \param name
  ///   Name of the configuration key whose value determines the concrete subclass of `PARAMETERS`
  ///   created during deserialization.
  /// \param description
  ///   Long description of this parameter.
  /// \param parent
  ///   Pointer to the Parameters object representing the collection of options located at
  ///   the same level of the configuration tree as `name`. A call to deserialize() or serialize()
  ///   on that object will automatically trigger a call to deserialize() or serialize() on this
  ///   parameter.
  explicit OptionalPolymorphicParameter(const char *name, const char * description,
                                        Parameters *parent)
    : ParameterBase(parent), name_(name), description_(description)
  {}

  OptionalPolymorphicParameter(const OptionalPolymorphicParameter &other)
    : name_(other.name_), description_(other.description_), id_(other.id_),
      value_(other.value_ ? other.value_->clone() : nullptr)
  {}

  OptionalPolymorphicParameter(OptionalPolymorphicParameter &&other) = default;

  OptionalPolymorphicParameter& operator=(const OptionalPolymorphicParameter &other) {
    OptionalPolymorphicParameter tmp(other);
    std::swap(*this, tmp);
    return *this;
  }

  OptionalPolymorphicParameter& operator=(OptionalPolymorphicParameter &&other) = default;

  ~OptionalPolymorphicParameter() override = default;

  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  void serialize(eckit::LocalConfiguration &config) const override;

  ObjectJsonSchema jsonSchema() const override;

  /// \brief Returns the identifier of the subclass of `PARAMETERS` whose instance is stored in this
  /// object, or boost::none if no value is stored.
  const boost::optional<std::string> &id() const { return id_; }

  /// \brief Pointer to the value stored in this parameter, or nullptr if no value is stored.
  const PARAMETERS *value() const { return value_.get(); }

  /// \brief Pointer to the value stored in this parameter, or nullptr if no value is stored.
  operator const PARAMETERS *() const { return value(); }

 private:
  typedef PolymorphicParameterTraits<PARAMETERS, FACTORY> Traits;

  std::string name_;
  std::string description_;
  boost::optional<std::string> id_;
  std::unique_ptr<PARAMETERS> value_;
};

template <typename PARAMETERS, typename FACTORY>
void OptionalPolymorphicParameter<PARAMETERS, FACTORY>::deserialize(
    util::CompositePath &path, const eckit::Configuration &config) {
  boost::optional<std::string> newId;
  std::unique_ptr<PARAMETERS> newValue;
  std::tie(newId, newValue) = Traits::get(path, config, name_);
  if (newId != boost::none && newValue != nullptr) {
    id_ = std::move(newId);
    value_ = std::move(newValue);
  }
}

template <typename PARAMETERS, typename FACTORY>
void OptionalPolymorphicParameter<PARAMETERS, FACTORY>::serialize(
    eckit::LocalConfiguration &config) const {
  if (id_ != boost::none && value_ != nullptr)
    Traits::set(config, name_, *id_, *value_);
}

template <typename PARAMETERS, typename FACTORY>
ObjectJsonSchema OptionalPolymorphicParameter<PARAMETERS, FACTORY>::jsonSchema() const {
  ObjectJsonSchema schema = Traits::jsonSchema(name_);
  if (description_ != "") {
    schema.extendPropertySchema(name_, {{"description", "\"" + description_ + "\""}});
  }
  return schema;
}

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_OPTIONALPOLYMORPHICPARAMETER_H_
