/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERTRAITS_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERTRAITS_H_

#include <string>

#include <boost/optional.hpp>

#include "eckit/config/Configuration.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

template <typename T>
struct ParameterTraits {
  /// \brief Retrieve the value of an option \p name from a configuration \p config.
  ///
  /// If the configuration contains an option with the required name and the conversion of its value
  /// into type T succeeds, the function returns that value. If the conversion fails, the function
  /// throws an exception, typically eckit::BadParameter. If the configuration does not contain
  /// an option with the required name, the function returns boost::none.
  static boost::optional<T> get(const eckit::Configuration &config,
                                const std::string& name) {
    T value;
    if (config.get(name, value))
      return value;
    else
      return boost::none;
  }
};

template <>
struct ParameterTraits<util::DateTime> {
  static boost::optional<util::DateTime> get(const eckit::Configuration &config,
                                             const std::string& name) {
    std::string value;
    if (config.get(name, value))
      return util::DateTime(value);
    else
      return boost::none;
  }
};

template <>
struct ParameterTraits<util::Duration> {
  static boost::optional<util::Duration> get(const eckit::Configuration &config,
                                             const std::string& name) {
    std::string value;
    if (config.get(name, value))
      return util::Duration(value);
    else
      return boost::none;
  }
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERTRAITS_H_
