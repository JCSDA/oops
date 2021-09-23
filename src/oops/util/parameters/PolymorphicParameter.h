/*
 * (C) Copyright 2020 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_POLYMORPHICPARAMETER_H_
#define OOPS_UTIL_PARAMETERS_POLYMORPHICPARAMETER_H_

#include <memory>
#include <string>
#include <tuple>
#include <utility>

#include <boost/optional.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ObjectJsonSchema.h"
#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/PolymorphicParameterTraits.h"

namespace oops {

/// \brief A parameter storing an instance of a subclass of `PARAMETERS` (itself a subclass of
/// Parameters).
///
/// The subclass is determined at runtime by loading the value of a specified key from a
/// Configuration object. If that key is not present, a default subclass is selected.
///
/// See RequiredPolymorphicParameter for general information about polymorphic parameters and an
/// explanation of the `FACTORY` template parameter.
///
/// \see RequiredPolymorphicParameter, OptionalPolymorphicParameter
template <typename PARAMETERS, typename FACTORY>
class PolymorphicParameter : public ParameterBase {
 public:
  /// \brief Constructor.
  ///
  /// \param name
  ///   Name of the configuration key whose value determines the concrete subclass of `PARAMETERS`
  ///   created during deserialization.
  /// \param defaultId
  ///   Identifier of the concrete subclass of `PARAMETERS` to be created if the above key is not
  ///   present.
  /// \param parent
  ///   Pointer to the Parameters object representing the collection of options located at
  ///   the same level of the configuration tree as `name`. A call to deserialize() on that object
  ///   will automatically trigger a call to deserialize() on the newly created object.
  PolymorphicParameter(const char *name, const char* defaultId, Parameters *parent)
    : PolymorphicParameter(name, "", defaultId, parent)
  {}

  /// \brief Constructor.
  ///
  /// \param name
  ///   Name of the configuration key whose value determines the concrete subclass of `PARAMETERS`
  ///   created during deserialization.
  /// \param description
  ///   Long description of this parameter.
  /// \param defaultId
  ///   Identifier of the concrete subclass of `PARAMETERS` to be created if the above key is not
  ///   present.
  /// \param parent
  ///   Pointer to the Parameters object representing the collection of options located at
  ///   the same level of the configuration tree as `name`. A call to deserialize() on that object
  ///   will automatically trigger a call to deserialize() on the newly created object.
  PolymorphicParameter(const char *name, const char *description, const char* defaultId,
                       Parameters *parent);

  PolymorphicParameter(const PolymorphicParameter &other)
    : name_(other.name_), description_(other.description_), id_(other.id_),
      value_(other.value_ ? other.value_->clone() : nullptr)
  {}

  PolymorphicParameter(PolymorphicParameter &&other) = default;

  PolymorphicParameter& operator=(const PolymorphicParameter &other) {
    PolymorphicParameter tmp(other);
    std::swap(*this, tmp);
    return *this;
  }

  PolymorphicParameter& operator=(PolymorphicParameter &&other) = default;

  ~PolymorphicParameter() override = default;

  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  void serialize(eckit::LocalConfiguration &config) const override;

  ObjectJsonSchema jsonSchema() const override;

  /// \brief Identifier of the subclass of `PARAMETERS` whose instance is stored in this object.
  const std::string &id() const { return id_; }

  /// \brief Value stored in this parameter.
  const PARAMETERS &value() const { return *value_; }

  /// \brief Value stored in this parameter.
  operator const PARAMETERS &() const { return value(); }

 private:
  typedef PolymorphicParameterTraits<PARAMETERS, FACTORY> Traits;

  std::string name_;
  std::string description_;
  std::string id_;
  std::string defaultId_;
  std::unique_ptr<PARAMETERS> value_;
};

template <typename PARAMETERS, typename FACTORY>
PolymorphicParameter<PARAMETERS, FACTORY>::PolymorphicParameter(const char *name,
                                                                const char *description,
                                                                const char* defaultId,
                                                                Parameters *parent)
  : ParameterBase(parent), name_(name), description_(description), defaultId_(defaultId) {
  util::CompositePath path;
  eckit::LocalConfiguration config;
  config.set(name, defaultId);
  deserialize(path, config);
}

template <typename PARAMETERS, typename FACTORY>
void PolymorphicParameter<PARAMETERS, FACTORY>::deserialize(util::CompositePath &path,
                                                            const eckit::Configuration &config) {
  boost::optional<std::string> newId;
  std::unique_ptr<PARAMETERS> newValue;
  std::tie(newId, newValue) = Traits::get(path, config, name_);
  if (newId != boost::none && newValue != nullptr) {
    id_ = std::move(*newId);
    value_ = std::move(newValue);
  }
}

template <typename PARAMETERS, typename FACTORY>
void PolymorphicParameter<PARAMETERS, FACTORY>::serialize(eckit::LocalConfiguration &config) const {
  Traits::set(config, name_, id_, *value_);
}

template <typename PARAMETERS, typename FACTORY>
ObjectJsonSchema PolymorphicParameter<PARAMETERS, FACTORY>::jsonSchema() const {
  ObjectJsonSchema schema = Traits::jsonSchema(name_);
  if (description_ != "") {
    schema.extendPropertySchema(name_, {{"description", "\"" + description_ + "\""}});
  }
  schema.extendPropertySchema(name_, {{"default", "\"" + defaultId_ + "\""}});
  return schema;
}

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_POLYMORPHICPARAMETER_H_
