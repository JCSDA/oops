/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERTRAITS_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERTRAITS_H_

#include <map>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/CompositePath.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/NamedEnumerator.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Traits dictating how parameters of type T are extracted from eckit::Configuration
/// objects.
///
/// Each specialization of this template must provide the following member functions:
///
/// \code
/// /// \brief Retrieve the value of the option \p name from the configuration \p config.
/// ///
/// /// If the configuration contains an option with the required name and the conversion of its
/// /// value into type T succeeds, the function returns that value. If the conversion fails, the
/// /// function throws an exception, typically eckit::BadParameter. If the configuration does not
/// /// contain an option with the required name, the function returns boost::none.
/// ///
/// /// \p path is the location of \p config in the full configuration loaded from a YAML file.
/// /// This object may be modified in-place by get(), but must be restored to its original state
/// /// before this function returns.
/// static boost::optional<T> get(util::CompositePath &path, const eckit::Configuration &config,
///                               const std::string &name);
/// \endcode
template <typename T,
          typename IsTDerivedFromParameters =
          typename std::integral_constant<bool,
                                          std::is_convertible<T*, Parameters*>::value>::type>
struct ParameterTraits;

/// \brief Generic implementation of the get() method for types for which this method
/// is overloaded in eckit::LocalConfiguration.
template <typename T>
struct GenericParameterTraits
{
  static boost::optional<T> get(util::CompositePath &path,
                                const eckit::Configuration &config,
                                const std::string &name) {
    T value;
    if (config.get(name, value)) {
      return value;
    } else {
      return boost::none;
    }
  }
};

/// \brief Can be used to define traits for parameters storing enums.
///
/// Example of use:
///
/// \code
/// enum class Fruit { APPLE, ORANGE };
///
/// struct FruitParameterTraitsHelper {
///   typedef Fruit EnumType;
///   static constexpr char enumTypeName[] = "Fruit";
///   static constexpr util::NamedEnumerator<Fruit> namedValues[] = {
///     { Fruit::APPLE, "apple" },
///     { Fruit::ORANGE, "orange" }
///   };
/// };
///
/// namespace oops {
///   template <>
///   struct ParameterTraits<Fruit> : public EnumParameterTraits<FruitParameterTraitsHelper> {};
/// }
///
/// // These two lines must be in a cc file. (They shouldn't be needed in C++17, which introduces
/// // inline variables.)
/// constexpr char FruitParameterTraitsHelper::enumTypeName[];
/// constexpr util::NamedEnumerator<Fruit> FruitParameterTraitsHelper::namedValues[];
/// \endcode
template <typename Helper>
struct EnumParameterTraits {
  typedef typename Helper::EnumType EnumType;
  static boost::optional<EnumType> get(util::CompositePath &path,
                                       const eckit::Configuration &config,
                                       const std::string &name) {
    std::string value;
    if (config.get(name, value)) {
      for (util::NamedEnumerator<EnumType> namedValue : Helper::namedValues)
        if (value == namedValue.name)
          return namedValue.value;

      util::PathComponent component(path, name);
      std::stringstream message;
      message << path.path() << ": The string '" << value << "' cannot be converted to "
              << Helper::enumTypeName;
      throw eckit::BadParameter(message.str(), Here());
    } else {
      return boost::none;
    }
  }
};

template <typename T>
struct IntegerParameterTraits : public GenericParameterTraits<T> {};

template <>
struct ParameterTraits<int, std::false_type> : public IntegerParameterTraits<int> {};

template <>
struct ParameterTraits<size_t, std::false_type> : public IntegerParameterTraits<size_t> {};

template <typename T>
struct FloatingPointParameterTraits : public GenericParameterTraits<T> {};

template <>
struct ParameterTraits<float, std::false_type> : public FloatingPointParameterTraits<float> {};

template <>
struct ParameterTraits<double, std::false_type> : public FloatingPointParameterTraits<double> {};

template <>
struct ParameterTraits<bool, std::false_type> : public GenericParameterTraits<bool> {};

template <>
struct ParameterTraits<std::string, std::false_type> :
    public GenericParameterTraits<std::string> {};

template <>
struct ParameterTraits<eckit::LocalConfiguration, std::false_type> :
    public GenericParameterTraits<eckit::LocalConfiguration> {};

/// \brief Specialization for types derived from Parameters.
template <typename T>
struct ParameterTraits<T, std::true_type>
{
  static boost::optional<T> get(util::CompositePath &path,
                                const eckit::Configuration &config,
                                const std::string& name) {
    if (config.has(name)) {
      util::PathComponent component(path, name);

      T value;
      value.deserialize(path, config.getSubConfiguration(name));
      return value;
    } else {
      return boost::none;
    }
  }
};

/// \brief Specialization for DateTime objects.
template <>
struct ParameterTraits<util::DateTime> {
  static boost::optional<util::DateTime> get(util::CompositePath &path,
                                             const eckit::Configuration &config,
                                             const std::string &name) {
    std::string value;
    if (config.get(name, value)) {
      return util::DateTime(value);
    } else {
      return boost::none;
    }
  }
};

/// \brief Specialization for Duration objects.
template <>
struct ParameterTraits<util::Duration> {
  static boost::optional<util::Duration> get(util::CompositePath &path,
                                             const eckit::Configuration &config,
                                             const std::string& name) {
    std::string value;
    if (config.get(name, value)) {
      return util::Duration(value);
    } else {
      return boost::none;
    }
  }
};

/// \brief Specialization for vectors.
template <typename Value>
struct ParameterTraits<std::vector<Value>, std::false_type> {
  static boost::optional<std::vector<Value>> get(util::CompositePath &path,
                                                 const eckit::Configuration &config,
                                                 const std::string &name) {
    boost::optional<std::vector<Value>> result;

    if (!config.has(name))
      return result;

    result = std::vector<Value>();

    std::string separator{config.separator()};
    for (const eckit::LocalConfiguration &subConfig : config.getSubConfigurations(name)) {
      util::PathComponent component(path, name);
      if (boost::optional<Value> value =
          ParameterTraits<Value>::get(path, subConfig, separator)) {
        result->push_back(std::move(*value));
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
  static boost::optional<std::map<Key, Value>> get(util::CompositePath &path,
                                                   const eckit::Configuration &config,
                                                   const std::string& name) {
    boost::optional<std::map<Key, Value>> result;

    if (!config.has(name))
      return result;

    result = std::map<Key, Value>();

    eckit::LocalConfiguration localConfig(config, name);
    for (const std::string &keyString : localConfig.keys()) {
      try {
        util::PathComponent component(path, keyString);
        Key key = boost::lexical_cast<Key>(keyString);
        if (boost::optional<Value> value =
            ParameterTraits<Value>::get(path, localConfig, keyString)) {
          (*result)[key] = *value;
        } else {
          throw eckit::Exception(path.path() + ": No value found for key '" + keyString +
                                 "'. Please note that because of a bug in the eckit YAML parser "
                                 "maps need to be written in the JSON style, with keys quoted. "
                                 "Example: my_int_to_float_map: {\"1\": 123, \"2\": 321}", Here());
        }
      }
      catch (boost::bad_lexical_cast &) {
        throw eckit::Exception(path.path() + ": Cannot convert '" + keyString + "' to key type",
                               Here());
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
