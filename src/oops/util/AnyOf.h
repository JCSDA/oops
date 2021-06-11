/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_ANYOF_H_
#define OOPS_UTIL_ANYOF_H_

#include <string>
#include <utility>

#include <boost/optional.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ObjectJsonSchema.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/TypeTraits.h"

namespace util {

/// \brief An object encapsulating the value of a Parameter whose type is unspecified until
/// runtime but must be one of the types listed in the `Types` template parameter.
///
/// Example: Suppose that a (mandatory) YAML configuration option `values` should be set to a list
/// of integers if certain conditions are met and to a string otherwise. You could then encapsulate
/// it in the following member variable of a Parameters subclass:
///
///     oops::RequiredParameter<util::AnyOf<std::string, std::vector<int>>> values{"values", this};
///
/// and retrieve the value of this option by calling the `as()` method:
///
///     if (/*some condition*/)
///       std::vector<int> v = values.value().as<std::vector<int>>();
///     else
///       std::string v = values.value().as<std::string>();
///
/// If the value of the YAML option is not in fact an instance of any of the types listed in the
/// `Types` parameter, validation of the YAML document against the JSON schema produced by the
/// Parameters object will fail. If the value of the YAML option is not an instance of the type
/// indicated in the call to `as()`, an exception will be thrown; the exception message will include
/// the path to the offending YAML node.
///
/// Note: to declare a `(Required|Optional)Parameter `object holding an `AnyOf<...>` value, include
/// not only the `AnyOf.h` header, but also the `ParameterTraitsAnyOf.h` header.
template <typename... Types>
class AnyOf {
 public:
  /// \brief Construct an AnyOf object holding the value \p initialValue.
  template <typename T>
  explicit AnyOf(const T &initialValue) : name_("dummy") {
    static_assert(util::any_is_same<T, Types...>::value,
                  "T must be one of the types from the Types list");
    oops::ParameterTraits<T>::set(config_, name_, initialValue);
  }

  /// \brief Construct an AnyOf object holding the value of the key \p name from the
  /// configuration \p config.
  ///
  /// \p path should be the location of \p config in the full configuration loaded from a YAML file.
  AnyOf(const util::CompositePath &path,
        const eckit::Configuration &config, const std::string &name)
    : path_(path), config_(config), name_(name)
  {
    ASSERT(config_.has(name_));
  }

  /// \brief Cast the stored value to type T and return it.
  ///
  /// An exception is thrown if the cast fails.
  template <typename T>
  T as() const {
    static_assert(util::any_is_same<T, Types...>::value,
                  "T must be one of the types from the Types list");
    util::CompositePath localPath = path_;
    boost::optional<T> result;
    try {
      // get() is expected to throw an exception if config_ has an option name_, but it cannot be
      // converted to type T.
      result = oops::ParameterTraits<T>::get(localPath, config_, name_);
    } catch (eckit::Exception &) {
      util::PathComponent component(localPath, name_);
      throw eckit::BadValue(localPath.path() + ": unexpected value type");
    }
    ASSERT(result != boost::none);
    return std::move(*result);
  }

  /// \brief Set the \p name key in the configuration \p config to the stored value.
  void serialize(eckit::LocalConfiguration &config, const std::string &name) const {
    eckit::LocalConfiguration value;
    config_.get(name_, value);  // read the value from the source Configuration...
    config.set(name, value);    // and store it in the destination Configuration.
  }

 private:
  util::CompositePath path_;
  eckit::LocalConfiguration config_;
  std::string name_;
};

}  // namespace util

#endif  // OOPS_UTIL_ANYOF_H_
