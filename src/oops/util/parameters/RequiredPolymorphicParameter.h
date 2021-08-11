/*
 * (C) Copyright 2020 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_REQUIREDPOLYMORPHICPARAMETER_H_
#define OOPS_UTIL_PARAMETERS_REQUIREDPOLYMORPHICPARAMETER_H_

#include <memory>
#include <string>
#include <tuple>
#include <utility>

#include <boost/optional.hpp>

#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/ObjectJsonSchema.h"
#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/PolymorphicParameterTraits.h"

namespace util {
class CompositePath;
}

namespace oops {

/// \brief A mandatory parameter storing an instance of a subclass of `PARAMETERS` (itself a
/// subclass of Parameters).
///
/// The subclass is determined at runtime by loading the value of a specified key from a
/// Configuration object. If that key is not present, an exception is thrown.
///
/// Polymorphic parameters represent branches of the configuration tree whose structure depends on
/// the value of a particular keyword. For example, here is a YAML file listing the properties of
/// some computer peripherals:
///
///     peripherals:
///       - type: mouse
///         has wheel: true
///       - type: printer
///         max page width (mm): 240
///         max page height (mm): 320
///
/// Clearly, the list of options that make sense for each item in the `peripherals` list depends on
/// the value of the `type` keyword. The structure of this YAML file could be represented with the
/// following subclasses of `Parameters`:
///
///     class PeripheralTypeDependentParameters : public Parameters {
///       OOPS_ABSTRACT_PARAMETERS(PeripheralTypeDependentParameters, Parameters)
///      public:
///       RequiredParameter<std::string> type{"type", this};
///     };
///
///     class PrinterParameters : public PeripheralTypeDependentParameters {
///       OOPS_CONCRETE_PARAMETERS(PrinterParameters, PeripheralTypeDependentParameters)
///      public:
///       RequiredParameter<int> maxPageWidth{"max page width", this};
///       RequiredParameter<int> maxPageHeight{"max page height", this};
///     };
///
///     class MouseParameters : public PeripheralTypeDependentParameters {
///       OOPS_CONCRETE_PARAMETERS(MouseParameters, PeripheralTypeDependentParameters)
///      public:
///       Parameter<int> numButtons{"num buttons", 3, this};
///     };
///
///     class PeripheralParameters : public Parameters {
///       OOPS_CONCRETE_PARAMETERS(PeripheralTypeDependentParameters, Parameters)
///      public:
///       RequiredPolymorphicParameter<PeripheralTypeDependentParameters, PeripheralFactory>
///         typeDependent{"type", this};
///     };
///
///     class ComputerParameters : public Parameters {
///       OOPS_CONCRETE_PARAMETERS(ComputerParameters, Parameters)
///      public:
///       Parameter<std::vector<PeripheralParameters>> peripherals{
///         "peripherals", {}, this};
///     };
///
/// Each item in the `peripherals` vector is represented with a
/// `RequiredPolymorphicParameter<PeripheralTypeDependentParameters, PeripheralFactory>` object.
/// This object holds a pointer to an instance of a subclass of the
/// `PeripheralTypeDependentParameters` abstract base class; whether it is an instance of
/// `PrinterParameters` or `MouseParameters` is determined at runtime depending on the value of the
/// `type` key. This is done by the `PeripheralFactory::createParameters()` static function (not
/// shown in the above code snippet), which takes the string loaded from the `type` key and returns
/// a unique pointer to a new instance of the subclass of `PeripheralTypeDependentParameters`
/// identified by that string. The `PeripheralFactory` class would typically be used also to create
/// objects representing the peripherals themselves.
///
/// \see PolymorphicParameter, OptionalPolymorphicParameter
template <typename PARAMETERS, typename FACTORY>
class RequiredPolymorphicParameter : public ParameterBase {
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
  explicit RequiredPolymorphicParameter(const char *name, Parameters *parent)
    : RequiredPolymorphicParameter(name, "", parent)
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
  explicit RequiredPolymorphicParameter(const char *name, const char *description,
                                        Parameters *parent)
    : ParameterBase(parent), name_(name), description_(description)
  {}

  RequiredPolymorphicParameter(const RequiredPolymorphicParameter &other)
    : name_(other.name_), description_(other.description_), id_(other.id_),
      value_(other.value_ ? other.value_->clone() : nullptr)
  {}

  RequiredPolymorphicParameter(RequiredPolymorphicParameter &&other) = default;

  RequiredPolymorphicParameter& operator=(const RequiredPolymorphicParameter &other) {
    RequiredPolymorphicParameter tmp(other);
    std::swap(*this, tmp);
    return *this;
  }

  RequiredPolymorphicParameter& operator=(RequiredPolymorphicParameter &&other) = default;

  ~RequiredPolymorphicParameter() override = default;

  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  void serialize(eckit::LocalConfiguration &config) const override;

  ObjectJsonSchema jsonSchema() const override;

  /// \brief Returns the identifier of the subclass of `PARAMETERS` whose instance is stored in this
  /// object.
  ///
  /// An exception is thrown if the value hasn't been loaded from a Configuration yet.
  const std::string &id() const {
    if (value_ == nullptr)
      throw std::runtime_error("The parameter hasn't been deserialized yet");
    return id_;
  }

  /// \brief Value stored in this parameter.
  ///
  /// An exception is thrown if the value hasn't been loaded from a Configuration yet.
  const PARAMETERS &value() const {
    if (value_ == nullptr)
      throw std::runtime_error("The parameter hasn't been deserialized yet");
    return *value_;
  }

  /// \brief Value stored in this parameter.
  ///
  /// An exception is thrown if the value hasn't been loaded from a Configuration yet.
  operator const PARAMETERS &() const { return value(); }

 private:
  typedef PolymorphicParameterTraits<PARAMETERS, FACTORY> Traits;

  std::string name_;
  std::string description_;
  std::string id_;
  std::unique_ptr<PARAMETERS> value_;
};

template <typename PARAMETERS, typename FACTORY>
void RequiredPolymorphicParameter<PARAMETERS, FACTORY>::deserialize(
    util::CompositePath &path, const eckit::Configuration &config) {
  boost::optional<std::string> newId;
  std::unique_ptr<PARAMETERS> newValue;
  std::tie(newId, newValue) = Traits::get(path, config, name_);
  if (newId == boost::none || newValue == nullptr)
    throw eckit::BadParameter(path.path() + ": Mandatory parameter '" + name_ + "' not found",
                              Here());
  id_ = std::move(*newId);
  value_ = std::move(newValue);
}

template <typename PARAMETERS, typename FACTORY>
void RequiredPolymorphicParameter<PARAMETERS, FACTORY>::serialize(
    eckit::LocalConfiguration &config) const {
  if (value_ != nullptr)
    Traits::set(config, name_, id_, *value_);
}

template <typename PARAMETERS, typename FACTORY>
ObjectJsonSchema RequiredPolymorphicParameter<PARAMETERS, FACTORY>::jsonSchema() const {
  ObjectJsonSchema schema = Traits::jsonSchema(name_);
  if (description_ != "") {
    schema.extendPropertySchema(name_, {{"description", "\"" + description_ + "\""}});
  }
  schema.require(name_);
  return schema;
}

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_REQUIREDPOLYMORPHICPARAMETER_H_
