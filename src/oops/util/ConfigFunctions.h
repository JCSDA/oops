/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_CONFIGFUNCTIONS_H_
#define OOPS_UTIL_CONFIGFUNCTIONS_H_

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

// ConfigFunctions are simple functions to manipulate eckit Configurations. These should be
// simple operations not tied to OOPS/JEDI concepts.

namespace util {

  /// \brief vectoriseAndFilter takes a single configuration object that is
  ///        internally a list of subConfigurations and breaks them up
  ///        into a vector of subConfigurations.
  ///        Then it loops over all subConfigurations and searches
  ///        for sub(sub)Configurations that are identified by the string \p tag.
  ///        If they exist, they are copied into the correct index of a
  ///        vector of LocalConfigurations.
  ///        If they don't exist an empty configuration is copied in the vector element.
  std::vector<eckit::LocalConfiguration>
    vectoriseAndFilter(const eckit::Configuration & config, const std::string & tag);

  /// \brief Check if the configuration is a vector of configurations
  bool isVector(const eckit::Configuration & config);

  /// \brief Check if the configuration is a subconfiguration
  bool isSubConfig(const eckit::Configuration & config);

  /// \brief Check if the configuration is a final pair
  bool isFinal(const eckit::Configuration & config);

  /// \brief Seek and replace a pattern with a value, recursively
  void seekAndReplace(eckit::LocalConfiguration & config, const std::string & pattern,
                      const std::string & value);
  void seekAndReplace(eckit::LocalConfiguration & config, const std::string & pattern,
                      const size_t & count, const size_t & zpad);

  /// \brief Merge configurations
  eckit::LocalConfiguration mergeConfigs(const eckit::Configuration & config1,
                                         const eckit::Configuration & config2);
}  // namespace util

#endif  // OOPS_UTIL_CONFIGFUNCTIONS_H_
