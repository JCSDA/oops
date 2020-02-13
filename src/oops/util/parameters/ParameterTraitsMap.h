/*
 * (C) Copyright 2020 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERTRAITSMAP_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERTRAITSMAP_H_

#include <map>
#include <string>

#include <boost/lexical_cast.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/ParameterTraits.h"

namespace oops {

/// \brief Traits needed for extraction of parameters storing maps.
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
struct ParameterTraits<std::map<Key, Value>> {
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
}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERTRAITSMAP_H_
