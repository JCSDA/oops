/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_OPTIONALPARAMETER_H_
#define OOPS_UTIL_PARAMETERS_OPTIONALPARAMETER_H_

#include <string>

#include <boost/optional.hpp>

#include "eckit/config/Configuration.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/ParameterTraits.h"

namespace oops {

/// \brief A parameter that does not need to be present in a Configuration and for which no
/// sensible default value can be defined.
template <typename T>
class OptionalParameter : public ParameterBase {
 public:
  explicit OptionalParameter(const char *name, Parameters *parent = nullptr)
    : ParameterBase(parent), name_(name)
  {}

  void deserialize(const eckit::Configuration &config) override {
    boost::optional<T> newValue = ParameterTraits<T>::get(config, name_);
    if (newValue != boost::none) {
      value_ = newValue;
    }
  }

  boost::optional<T> value() const { return value_; }

  operator boost::optional<T>() const { return value_; }

 private:
  std::string name_;
  boost::optional<T> value_;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_OPTIONALPARAMETER_H_
