
/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/ConfigFunctions.h"

namespace util {

  namespace oopsconfigfunctions {

    std::vector<eckit::LocalConfiguration>
       vectoriseAndFilter(const eckit::Configuration & config, const std::string & tag)
    {
      const std::vector<eckit::LocalConfiguration>
        ObsConfigVec(config.getSubConfigurations());
      std::vector<eckit::LocalConfiguration> filteredConfig;
      for (auto conf : ObsConfigVec) {
          eckit::LocalConfiguration tempConf;
          if (conf.has(tag)) {
            tempConf = eckit::LocalConfiguration(conf, tag);
          }
          filteredConfig.push_back(tempConf);
      }
      return filteredConfig;
    }

  }  // namespace oopsconfigfunctions

}  // namespace util
