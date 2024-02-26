/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERTRAITS_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERTRAITS_H_

#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
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
#include "oops/util/parameters/ObjectJsonSchema.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/PartialDateTime.h"
#include "oops/util/stringFunctions.h"  // for join()

namespace oops {

/// \brief Traits dictating how parameters of type T are transferred to and from
/// eckit::Configuration objects.
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
///
/// /// \brief Save the value \p value of the option \p name to the configuration \p config.
/// static void set(eckit::LocalConfiguration &config, const std::string &name, const T &value);
///
/// /// \brief Return an object encapsulating the JSON schema specifying the expected structure of
/// /// the YAML node from which option \p name is loaded.
/// static ObjectJsonSchema jsonSchema(const std::string &name);
/// \endcode
///
/// Optionally, the following member function may also be provided:
///
/// \code
/// /// \brief Return a JSON representation of a value of type \c T.
/// static std::string valueAsJson(const T &value);
/// \endcode
///
/// If this member function is available, default values of parameters holding values of type \c T
/// can be embedded in JSON schemas.
template <typename T,
          typename IsTDerivedFromParameters =
          typename std::integral_constant<bool,
                                          std::is_convertible<T*, Parameters*>::value>::type>
struct ParameterTraits;

/// \brief Generic implementation of the get() and set() methods for types for which these methods
/// are overloaded in eckit::LocalConfiguration.
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

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const T &value) {
    config.set(name, value);
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

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const EnumType &value) {
    for (util::NamedEnumerator<EnumType> namedValue : Helper::namedValues)
      if (value == namedValue.value) {
        config.set(name, namedValue.name);
        return;
      }

    std::stringstream message;
    message << "The value '" << std::to_string(static_cast<std::intmax_t>(value)) <<
               "' cannot be converted to a " << Helper::enumTypeName;
    throw eckit::BadParameter(message.str(), Here());
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    std::string enumSchema = '[' +
        util::stringfunctions::join(", ",
                                    std::begin(Helper::namedValues),
                                    std::end(Helper::namedValues),
                                    quotedName) +
        ']';
    return ObjectJsonSchema({{name, {{"enum", enumSchema}}}});
  }

  static std::string valueAsJson(const EnumType &value) {
    for (util::NamedEnumerator<EnumType> namedValue : Helper::namedValues)
      if (value == namedValue.value) {
        return quotedName(namedValue);
      }
    return "";
  }

 private:
  static std::string quotedName(util::NamedEnumerator<EnumType> namedValue) {
    std::string result = "\"";
    result += namedValue.name;
    result += "\"";
    return result;
  }
};

template <typename T>
struct IntegerParameterTraits : public GenericParameterTraits<T> {
  static ObjectJsonSchema jsonSchema(const std::string &name) {
    return ObjectJsonSchema({{name, {{"type", "\"integer\""}}}});
  }

  static std::string valueAsJson(const T &value);
};

template <>
struct ParameterTraits<int, std::false_type> : public IntegerParameterTraits<int> {};

template <>
struct ParameterTraits<size_t, std::false_type> : public IntegerParameterTraits<size_t> {};

template <>
struct ParameterTraits<int64_t, std::false_type> : public IntegerParameterTraits<int64_t> {};

template <typename T>
struct FloatingPointParameterTraits : public GenericParameterTraits<T> {
  static ObjectJsonSchema jsonSchema(const std::string &name) {
    return ObjectJsonSchema({{name, {{"type", "\"number\""}}}});
  }

  static std::string valueAsJson(const T &value);
};

template <>
struct ParameterTraits<float, std::false_type> : public FloatingPointParameterTraits<float> {};

template <>
struct ParameterTraits<double, std::false_type> : public FloatingPointParameterTraits<double> {};

template <>
struct ParameterTraits<bool, std::false_type> : public GenericParameterTraits<bool> {
  static ObjectJsonSchema jsonSchema(const std::string &name) {
    return ObjectJsonSchema({{name, {{"type", "\"boolean\""}}}});
  }

  static std::string valueAsJson(const bool &value);
};

template <>
struct ParameterTraits<std::string, std::false_type> :
    public GenericParameterTraits<std::string> {
  static ObjectJsonSchema jsonSchema(const std::string &name) {
    return ObjectJsonSchema({{name, {{"type", R"(["string", "number"])"}}}});
  }

  static std::string valueAsJson(const std::string &value);
};

template <>
struct ParameterTraits<eckit::LocalConfiguration, std::false_type> :
    public GenericParameterTraits<eckit::LocalConfiguration> {
  static ObjectJsonSchema jsonSchema(const std::string &name) {
    return ObjectJsonSchema({{name, PropertyJsonSchema()}});
  }
};

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

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const T &value) {
    eckit::LocalConfiguration subConfig;
    value.serialize(subConfig);
    config.set(name, subConfig);
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    T value;
    return ObjectJsonSchema({{name, value.jsonSchema().toPropertyJsonSchema()}});
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

  static void set(eckit::LocalConfiguration &config,
             const std::string &name,
             const util::DateTime &value) {
    config.set(name, value.toString());
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    return ObjectJsonSchema({{name, {{"type", "\"string\""},
                                     {"format", "\"date-time\""}}}});
  }

  static std::string valueAsJson(const util::DateTime &value);
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

  static void set(eckit::LocalConfiguration &config,
             const std::string &name,
             const util::Duration &value) {
    config.set(name, value.toString());
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    return ObjectJsonSchema({{name, {{"type", "\"string\""},
                                     {"format", "\"duration\""}}}});
  }

  static std::string valueAsJson(const util::Duration &value);
};


/// \brief Specialization for PartialDateTime objects.
template <>
struct ParameterTraits<util::PartialDateTime> {
  static boost::optional<util::PartialDateTime> get(util::CompositePath &path,
                                                    const eckit::Configuration &config,
                                                    const std::string& name) {
    std::string value;
    if (config.get(name, value)) {
      return util::PartialDateTime(value);
    } else {
      return boost::none;
    }
  }

  static void set(eckit::LocalConfiguration &config,
             const std::string &name,
             const util::PartialDateTime &value) {
    config.set(name, value.toString());
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    std::string expression = R"("^([0-9]{4}|\\*{4})-([0-9]{2}|\\*{2})-([0-9]{2}|\\*{2}))"
                             R"(T([0-9]{2}|\\*{2}):([0-9]{2}|\\*{2}):([0-9]{2}|\\*{2})Z$")";
    return ObjectJsonSchema({{name, {{"type", "\"string\""},
                                     {"pattern", expression}}}});
  }

  static std::string valueAsJson(const util::PartialDateTime &value);
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

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const std::vector<Value> &value) {
    std::vector<eckit::LocalConfiguration> subConfigs;
    subConfigs.reserve(value.size());

    eckit::LocalConfiguration temp;
    const std::string dummyName = "a";
    for (size_t i = 0; i < value.size(); ++i) {
      ParameterTraits<Value>::set(temp, dummyName, value[i]);
      subConfigs.push_back(temp.getSubConfiguration(dummyName));
    }
    config.set(name, subConfigs);
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    ObjectJsonSchema itemSchema = ParameterTraits<Value>::jsonSchema("");
    return ObjectJsonSchema({{name, {{"type", "\"array\""},
                                     {"items", toString(itemSchema.properties().at(""))}}}});
  }
};

/// \brief Specialization for pairs (specified as a two-element list)
template <typename Value1, typename Value2>
struct ParameterTraits<std::pair<Value1, Value2>, std::false_type> {
  static boost::optional<std::pair<Value1, Value2>> get(util::CompositePath &path,
                                                        const eckit::Configuration &config,
                                                        const std::string &name) {
    boost::optional<std::pair<Value1, Value2>> result;

    if (!config.has(name))
      return result;

    const std::vector<eckit::LocalConfiguration> subConfigs = config.getSubConfigurations(name);
    if (subConfigs.size() != 2) {
      throw eckit::Exception(path.path() + " has to have exactly two entries for key '" +
                             name + "'", Here());
    }
    util::PathComponent component(path, name);
    std::string separator{config.separator()};
    boost::optional<Value1> value1 = ParameterTraits<Value1>::get(path, subConfigs[0], separator);
    boost::optional<Value2> value2 = ParameterTraits<Value2>::get(path, subConfigs[1], separator);
    if (value1 && value2) {
      result = std::pair<Value1, Value2>(std::move(*value1), std::move(*value2));
    }
    return result;
  }

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const std::pair<Value1, Value2> &values) {
    std::vector<eckit::LocalConfiguration> subConfigs;
    subConfigs.reserve(2);

    eckit::LocalConfiguration temp;
    const std::string dummyName = "a";
    ParameterTraits<Value1>::set(temp, dummyName, values.first);
    subConfigs.push_back(temp.getSubConfiguration(dummyName));
    ParameterTraits<Value2>::set(temp, dummyName, values.second);
    subConfigs.push_back(temp.getSubConfiguration(dummyName));

    config.set(name, subConfigs);
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    std::stringstream items;
    {
      eckit::Channel ch;
      ch.setStream(items);
      ch << "[\n";
      {
        eckit::AutoIndent indent(ch);
        ObjectJsonSchema item1Schema = ParameterTraits<Value1>::jsonSchema("");
        ch << toString(item1Schema.properties().at(""));
        ch << ",\n";
        ObjectJsonSchema item2Schema = ParameterTraits<Value2>::jsonSchema("");
        ch << toString(item2Schema.properties().at(""));
        ch << '\n';
      }
      ch << "]";
    }

    return ObjectJsonSchema({{name, {{"type", "\"array\""},
                                     {"items", items.str()},
                                     {"minItems", "2"},
                                     {"maxItems", "2"}}}});
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

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const std::map<Key, Value> &value) {
    for (const std::pair<const Key, Value> &keyValue : value)
      ParameterTraits<Value>::set(
        config,
        name + config.separator() + boost::lexical_cast<std::string>(keyValue.first),
        keyValue.second);
    // If the map is empty, the loop above won't do anything, so we need to set the 'name' key
    // to an empty value of type "map" explicitly.
    if (value.empty())
      config.set(name, eckit::LocalConfiguration());
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    std::stringstream patternProperties;
    {
      eckit::Channel ch;
      ch.setStream(patternProperties);
      ch << "{\n";
      {
        eckit::AutoIndent indent(ch);
        ch << R"("" : )";
        ObjectJsonSchema itemSchema = ParameterTraits<Value>::jsonSchema("");
        ch << toString(itemSchema.properties().at(""));
        ch << '\n';
      }
      ch << "}";
    }

    return ObjectJsonSchema({{name, {{"type", "\"object\""},
                                     {"patternProperties", patternProperties.str()}}}});
  }
};

/// \brief Specialization for set<int>.
///
/// This specialization handles conversion of std::set<int> to and from string-valued YAML options
/// having the form of comma-separated lists of single integers or ranges of integers. Examples:
///
///     option_a: 1
///     option_b: 10-15
///     option_c: 8,10, 11-15, 1,3
///
/// \warning The current implementation can only handle non-negative integers.
template <>
struct ParameterTraits<std::set<int>, std::false_type> {
  static boost::optional<std::set<int>> get(util::CompositePath &path,
                                            const eckit::Configuration &config,
                                            const std::string &name);

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const std::set<int> &value);

  static ObjectJsonSchema jsonSchema(const std::string &name);

  static std::string valueAsJson(const std::set<int> &value);
};


// Note: to avoid coupling this file too tightly with headers defining a lot of disparate classes,
// only specializations of ParameterTraits for commonly used types should be added here.
// Specializations for types used less frequently should be defined in separate files (see
// ParameterTraitsScalarOrMap.h for an example).

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERTRAITS_H_
