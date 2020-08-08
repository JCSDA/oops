/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERTRAITS_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERTRAITS_H_

#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/ParameterBase.h"

namespace oops {

/// \brief Traits dictating how parameters of type T are extracted from eckit::Configuration
/// objects.
template <typename T,
          typename IsTDerivedFromParameters =
          typename std::integral_constant<bool,
                                          std::is_convertible<T*, Parameters*>::value>::type>
struct ParameterTraits;

/// \brief Generic implementation for types not derived from Parameters.
template <typename T>
struct ParameterTraits<T, std::false_type>
{
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

/// \brief Specialization for types derived from Parameters.
template <typename T>
struct ParameterTraits<T, std::true_type>
{
  static boost::optional<T> get(const eckit::Configuration &config, const std::string& name) {
    if (config.has(name)) {
      T value;
      value.deserialize(config.getSubConfiguration(name));
      return value;
    } else {
      return boost::none;
    }
  }
};

/// \brief Specialization for DateTime objects.
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

/// \brief Specialization for Duration objects.
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

/// \brief Specialization for vectors.
template <typename Value>
struct ParameterTraits<std::vector<Value>, std::false_type> {
  static boost::optional<std::vector<Value>> get(const eckit::Configuration &config,
                                                 const std::string& name) {
    boost::optional<std::vector<Value>> result;

    if (!config.has(name))
      return result;

    result = std::vector<Value>();

    std::string separator{config.separator()};
    for (const eckit::LocalConfiguration &subConfig : config.getSubConfigurations(name)) {
      if (boost::optional<Value> value =
          ParameterTraits<Value>::get(subConfig, separator)) {
        result->push_back(*value);
      }
    }

    return result;
  }
};

/// \brief Specialization for maps.
///
/// \note Owing to a bug in the eckit YAML parser, maps need to be written in the JSON style,
/// with keys quoted. Example:
///
///   my_int_to_float_map: {"1": 123, "2": 321}
///
/// This will not work:
///
///   my_int_to_float_map: {1: 123, 2: 321}
///
/// and neither will this:
///
///   my_int_to_float_map:
///     1: 123
///     2: 321
///
/// or this:
///
///   my_int_to_float_map:
///     "1": 123
///     "2": 321
template <typename Key, typename Value>
struct ParameterTraits<std::map<Key, Value>, std::false_type> {
  static boost::optional<std::map<Key, Value>> get(const eckit::Configuration &config,
                                                   const std::string& name) {
    boost::optional<std::map<Key, Value>> result;

    if (!config.has(name))
      return result;

    result = std::map<Key, Value>();

    eckit::LocalConfiguration localConfig(config, name);
    for (const std::string &keyString : localConfig.keys()) {
      try {
        Key key = boost::lexical_cast<Key>(keyString);
        if (boost::optional<Value> value = ParameterTraits<Value>::get(localConfig, keyString)) {
          (*result)[key] = *value;
        } else {
          throw eckit::Exception("No value found for key '" + keyString +
                                 "'. Please note that because of a bug in the eckit YAML parser "
                                 "maps need to be written in the JSON style, with keys quoted. "
                                 "Example: my_int_to_float_map: {\"1\": 123, \"2\": 321}", Here());
        }
      }
      catch (boost::bad_lexical_cast &) {
        throw eckit::Exception("Cannot convert '" + keyString + "' to key type", Here());
      }
    }

    return result;
  }
};

// Note: to avoid coupling this file too tightly with headers defining a lot of disparate classes,
// only specializations of ParameterTraits for commonly used types should be added here.
// Specializations for types used less frequently should be defined in separate files (see
// ParameterTraitsScalarOrMap.h for an example).

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERTRAITS_H_
