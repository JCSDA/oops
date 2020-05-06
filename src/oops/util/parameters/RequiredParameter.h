/*
 * (C) Copyright 2020 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_REQUIREDPARAMETER_H_
#define OOPS_UTIL_PARAMETERS_REQUIREDPARAMETER_H_

#include <set>
#include <string>

#include <boost/optional.hpp>

#include "eckit/config/Configuration.h"
#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/ParameterTraits.h"

namespace oops {

/// \brief A parameter that must be present in a Configuration (otherwise an exception is thrown).
template <typename T>
class RequiredParameter : public ParameterBase {
 public:
  explicit RequiredParameter(const char *name, Parameters *parent = nullptr)
    : ParameterBase(parent), name_(name)
  {}

  void deserialize(const eckit::Configuration &config, std::set<std::string> &usedKeys) override {
    boost::optional<T> newValue = ParameterTraits<T>::get(config, name_);
    if (newValue == boost::none) {
      throw eckit::BadParameter("Mandatory parameter '" + name_ + "' not found", Here());
    }
    value_ = newValue;
    usedKeys.insert(name_);
  }

  /// \brief Returns the stored value.
  ///
  /// An exception is thrown if the value hasn't been loaded from a Configuration yet.
  const T &value() const { return value_.value(); }

  /// \brief Returns the stored value.
  ///
  /// An exception is thrown if the value hasn't been loaded from a Configuration yet.
  operator const T &() const { return value_.value(); }

 private:
  std::string name_;
  // The value is stored in a boost::optional object because T may not be
  // default-constructible.
  boost::optional<T> value_;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_REQUIREDPARAMETER_H_
