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
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Traits dictating how instances of subclasses of `PARAMETERS` (itself a subclass of
/// Parameters) are transferred to and from eckit::Configuration objects.
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
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_POLYMORPHICPARAMETERTRAITS_H_
