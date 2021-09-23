/*
 * (C) Copyright 2021 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERTRAITSANYOF_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERTRAITSANYOF_H_

#include <string>
#include <utility>
#include <vector>

#include <boost/optional.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/AnyOf.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/stringFunctions.h"  // for join()

namespace oops {

/// \brief Traits dictating how parameters of type AnyOf<...> are transferred to and from
/// eckit::Configuration objects.

template <typename... Types>
struct ParameterTraits<util::AnyOf<Types...>, std::false_type> {
  static boost::optional<util::AnyOf<Types...>> get(util::CompositePath &path,
                                                    const eckit::Configuration &config,
                                                    const std::string &name) {
    boost::optional<util::AnyOf<Types...>> result;

    if (!config.has(name))
      return result;

    result = util::AnyOf<Types...>(path, config, name);
    return result;
  }

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const util::AnyOf<Types...> &value) {
    value.serialize(config, name);
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    const std::vector<ObjectJsonSchema> typeSchemas{ParameterTraits<Types>::jsonSchema("")...};
    std::string anyOfSchema = '[' +
        util::stringfunctions::join(", ", typeSchemas.begin(), typeSchemas.end(),
                                    [](const ObjectJsonSchema &s)
                                    {return toString(s.properties().at(""));}) +
        ']';
    return ObjectJsonSchema({{name, {{"anyOf", std::move(anyOfSchema)}}}});
  }

  static ObjectJsonSchema jsonSchema(const std::string &name,
                                     const util::AnyOf<Types...> &defaultValue)
  {
    return jsonSchema(name);
  }
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERTRAITSANYOF_H_
