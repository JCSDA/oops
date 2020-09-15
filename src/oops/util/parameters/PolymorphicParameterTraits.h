/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_POLYMORPHICPARAMETERTRAITS_H_
#define OOPS_UTIL_PARAMETERS_POLYMORPHICPARAMETERTRAITS_H_

#include <map>
#include <memory>
#include <string>
#include <utility>

#include <boost/optional.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ObjectJsonSchema.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Traits dictating how instances of subclasses of `PARAMETERS` (itself a subclass of
/// Parameters) are transferred to and from eckit::Configuration objects.
///
/// \tparam PARAMETERS
///   Common base class of the classes whose transfer to and from Configuration objects is
///   controlled by this instantiation of PolymorphicParameterTraits.
///
/// \tparam FACTORY
///   A factory class able to create instances of subclasses of `PARAMETERS` and enumerate
///   valid identifiers of these subclasses. It needs to provide the following functions:
///
///       /// Return a new instance of the subclass of PARAMETERS identified by the string `id`.
///       static std::unique_ptr<PARAMETERS> createParameters(const std::string &id);
///
///       /// Return the list of identifiers recognized by createParameters().
///       ///
///       /// Here, StringContainer can be any container of `std::string` objects that can be
///       /// iterated over by a range-based `for` loop (e.g. `std::vector<std::string>` or
///       /// `std::set<std::string>`); it can be returned by value or by reference.
///       static StringContainer getMakerNames();
template <typename PARAMETERS, typename FACTORY>
struct PolymorphicParameterTraits {
  /// \brief Try to deserialize a polymorphic parameter from a Configuration object.
  ///
  /// \returns If \p config has a key \p name, a pair composed of the value of that key
  /// and a pointer to an instance of a subclass of PARAMETERS created by passing the
  /// value of that key to FACTORY::createParameters(). If key \p name is not present,
  /// the pair `(boost::none, nullptr)`.
  static std::pair<boost::optional<std::string>, std::unique_ptr<PARAMETERS>> get(
      util::CompositePath &path,
      const eckit::Configuration &config,
      const std::string& name) {
    std::pair<boost::optional<std::string>, std::unique_ptr<PARAMETERS>> result;
    if (config.has(name)) {
      result.first = config.getString(name);
      result.second = FACTORY::createParameters(*result.first);
      result.second->deserialize(path, config);
    }
    return result;
  }

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const std::string &id,
                  const PARAMETERS &value) {
    config.set(name, id);
    value.serialize(config);
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    std::map<std::string, ObjectJsonSchema> variants;
    for (const std::string &id : FACTORY::getMakerNames()) {
      std::unique_ptr<PARAMETERS> params = FACTORY::createParameters(id);
      variants[id] = params->jsonSchema();
    }
    return ObjectJsonSchema(name, variants);
  }
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_POLYMORPHICPARAMETERTRAITS_H_
