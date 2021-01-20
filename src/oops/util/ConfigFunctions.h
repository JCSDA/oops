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


namespace util {
  /// Namespace for oopsconfigfunctions
  namespace oopsconfigfunctions {

    /// \brief vectoriseAndFilter takes a single configuration object that is
    ///        internally a list of subConfigurations and breaks them up
    ///        into a vector of subConfigurations.
    ///        Then it loops over all subConfigurations and searches
    ///        for sub(sub)Configurations that are identified by the string tag.
    ///        If they exist, they are copied into the correct index of a
    ///        vector of LocalConfigurations.
    ///        If they don't exist an empty configuration is copied in the vector element.
    std::vector<eckit::LocalConfiguration>
      vectoriseAndFilter(const eckit::Configuration &config, const std::string & tag);

  }  // namespace oopsconfigfunctions
}  // namespace util

#endif  // OOPS_UTIL_CONFIGFUNCTIONS_H_
