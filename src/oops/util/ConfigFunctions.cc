
/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/ConfigFunctions.h"

namespace util {

  std::vector<eckit::LocalConfiguration>
    vectoriseAndFilter(const eckit::Configuration & config, const std::string & tag)
  {
    const std::vector<eckit::LocalConfiguration>
      ObsConfigVec(config.getSubConfigurations());
    std::vector<eckit::LocalConfiguration> filteredConfig;
    for (const eckit::LocalConfiguration & conf : ObsConfigVec) {
      eckit::LocalConfiguration tempConf = conf.getSubConfiguration(tag);
      filteredConfig.push_back(tempConf);
    }
    return filteredConfig;
  }

}  // namespace util
