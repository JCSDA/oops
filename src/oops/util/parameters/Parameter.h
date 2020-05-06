/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETER_H_
#define OOPS_UTIL_PARAMETERS_PARAMETER_H_

#include <set>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraits.h"

namespace oops {

/// \brief A parameter with a default value.
///
/// \note Before declaring a variable of type Parameter<T> where T is something else than a plain
/// old type (int, float etc.), std::map, std::string, std::vector, util::DateTime, util::Duration
/// or a class derived from ParameterBase, check if there is a dedicated ParameterTraits*.h header
/// containing the specialization of ParameterTraits for type T. If so, include it before the
/// variable declaration. (The "main" ParameterTraits.h header only contains specializations of
/// ParameterTraits for the types listed above).
template <typename T>
class Parameter : public ParameterBase {
 public:
  Parameter(const char *name, const T& defaultValue, Parameters *parent = nullptr)
    : ParameterBase(parent), name_(name), value_(defaultValue)
  {}

  void deserialize(const eckit::Configuration &config, std::set<std::string> &usedKeys) override {
    boost::optional<T> newValue = ParameterTraits<T>::get(config, name_);
    if (newValue != boost::none) {
      value_ = newValue.get();
      usedKeys.insert(name_);
    }
  }

  const T &value() const { return value_; }

  operator const T &() const { return value_; }

 private:
  std::string name_;
  T value_;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETER_H_
