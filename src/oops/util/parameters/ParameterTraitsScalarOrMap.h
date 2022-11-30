/*
 * (C) Copyright 2020 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERTRAITSSCALARORMAP_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERTRAITSSCALARORMAP_H_

#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/utils/StringTools.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/ScalarOrMap.h"

namespace oops {

/// \brief Traits needed for extraction of parameters storing either maps or scalars.
///
/// \note Maps need to be written in the JSON style; see the documentation of
/// ParameterTraits<std::map<Key, Value>> for more details.
template <typename Key, typename Value>
struct ParameterTraits<util::ScalarOrMap<Key, Value>, std::false_type>
{
  static boost::optional<util::ScalarOrMap<Key, Value>> get(util::CompositePath &path,
                                                            const eckit::Configuration &config,
                                                            const std::string& name) {
    if (!config.has(name))
      return boost::none;

    std::list<std::string> messages;
    try {
      // Try to extract a scalar from config.name.
      if (boost::optional<Value> scalar =
          ParameterTraits<Value>::get(path, config, name)) {
        return util::ScalarOrMap<Key, Value>(*scalar);
      }
    }
    catch (eckit::Exception &e) {
      // It isn't a scalar.
      messages.push_back(e.what());
    }

    try {
      // Try to extract a map from config.name.
      if (boost::optional<std::map<Key, Value>> map =
          ParameterTraits<std::map<Key, Value>>::get(path, config, name)) {
        return util::ScalarOrMap<Key, Value>(*map);
      }
    }
    catch (eckit::Exception &e) {
      // It isn't a map.
      messages.push_back(e.what());
    }

    messages.push_front("The contents of '" + name + "' are neither a scalar nor a map.");
    throw eckit::Exception(eckit::StringTools::join("\n", messages.begin(), messages.end()),
                           Here());
  }

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const util::ScalarOrMap<Key, Value> &value) {
    if (value.isScalar()) {
      ParameterTraits<Value>::set(config, name, value.at(Key()));
    } else {
      for (const std::pair<const Key, Value> &keyValue : value)
        ParameterTraits<Value>::set(
          config,
          name + config.separator() + boost::lexical_cast<std::string>(keyValue.first),
          keyValue.second);
      // If the map is empty, the loop above won't do anything, so we need to set the 'name' key
      // to an empty value of type "map" explicitly.
      if (value.begin() == value.end())  // is the map empty?
        config.set(name, eckit::LocalConfiguration());
    }
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    std::stringstream oneOf;
    {
      eckit::Channel ch;
      ch.setStream(oneOf);
      ch << "[\n";
      {
        eckit::AutoIndent indent(ch);
        ObjectJsonSchema scalarSchema = ParameterTraits<Value>::jsonSchema("");
        ObjectJsonSchema mapSchema = ParameterTraits<std::map<Key, Value>>::jsonSchema("");
        ch << toString(scalarSchema.properties().at("")) << ",\n";
        ch << toString(mapSchema.properties().at("")) << '\n';
      }
      ch << "]";
    }

    return ObjectJsonSchema({{name, {{"oneOf", oneOf.str()}}}});
  }

  static ObjectJsonSchema jsonSchema(const std::string &name,
                                     const util::ScalarOrMap<Key, Value> &defaultValue)
  {
    return jsonSchema(name);
  }
};
}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERTRAITSSCALARORMAP_H_
