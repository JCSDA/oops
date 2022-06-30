
/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/ConfigFunctions.h"

#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>

#include "eckit/log/JSON.h"
#include "oops/util/abor1_cpp.h"

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

  bool isVector(const eckit::Configuration & config)
  {
    // Print configuration into stringstream, then into string
    std::stringstream ss;
    ss << config << std::endl;
    std::string str = ss.str();

    // Check if the configuration is a vector of configurations
    return str.find("LocalConfiguration[root=(") != std::string::npos;
  }

  bool isSubConfig(const eckit::Configuration & config)
  {
    // Print configuration into stringstream, then into string
    std::stringstream ss;
    ss << config << std::endl;
    std::string str = ss.str();

    // Check if the configuration is a subconfiguration
    return str.find("LocalConfiguration[root={") != std::string::npos;
  }

  bool isFinal(const eckit::Configuration & config)
  {
    // Check if the configuration is a final pair
    return !(isVector(config) || isSubConfig(config));
  }

  void seekAndReplace(eckit::LocalConfiguration & config, const std::string & pattern,
                      const std::string & value_out)
  {
    // Check if config is a subconfiguration
    if (isSubConfig(config)) {
      // Get the subconfiguration keys
      std::vector<std::string> keys(config.keys());

      // Loop over keys
      for (size_t jj = 0; jj < keys.size(); ++jj) {
        // Get value for the key "keys[jj]" as a local configuration
        eckit::LocalConfiguration subConfig(config, keys[jj]);

        // Check this local configuration nature
        if (isVector(subConfig)) {
          // The local configuration is a vector of configurations

          // Get items as local configurations
          std::vector<eckit::LocalConfiguration> subConfigs;
          config.get(keys[jj], subConfigs);

          // Check items nature
          if (isSubConfig(subConfigs[0])) {
            // Items are subconfigurations
            std::vector<eckit::LocalConfiguration> updatedSubConfigs;

            // Loop over items
            for (eckit::LocalConfiguration & subConfig : subConfigs) {
              // Call seekAndReplace for a subconfiguration
              seekAndReplace(subConfig, pattern, value_out);
              updatedSubConfigs.push_back(subConfig);
            }

            // Reset vector
            config.set(keys[jj], updatedSubConfigs);
          } else if (isVector(subConfigs[0])) {
            // Items are vectors of subconfigurations
            ABORT("Vector of vectors is not implemented yet...");
          } else {
            // Items are final pairs

            // Get items as strings
            std::vector<std::string> subStrings;
            config.get(keys[jj], subStrings);

            // Loop over items
            for (size_t kk=0; kk < subStrings.size(); ++kk) {
              // Check if the substring contains the pattern
              if (subStrings[kk].find(pattern) != std::string::npos) {
                // Update the substring
                subStrings[kk] = std::regex_replace(subStrings[kk], std::regex(pattern), value_out);
              }
            }

            // Reset vector of final pairs
            config.set(keys[jj], subStrings);
          }
        } else if (isFinal(subConfig)) {
          // The local configuration is a final pair

          // Get value for the key "keys[jj]" as a string
          std::string value(config.getString(keys[jj]));

          // Check if the value contains the pattern
          if (value.find(pattern) != std::string::npos) {
            // Update the value
            value = std::regex_replace(value, std::regex(pattern), value_out);

            // Reset the value
            config.set(keys[jj], value);
          }
        } else {
          // The local configuration is itself a subconfiguration

          // Call seekAndReplace for a subconfiguration
          seekAndReplace(subConfig, pattern, value_out);

          // Reset subconfiguration
          config.set(keys[jj], subConfig);
        }
      }
    } else if (isVector(config)) {
       // seekAndReplace should not be called on a vector of configurations
       ABORT("This is a vector, shouldn't be here...");
    } else {
       // seekAndReplace should not be called on a final pair
       ABORT("This is a final pair, shouldn't be here...");
    }
  }

  void seekAndReplace(eckit::LocalConfiguration & config, const std::string & pattern,
                      const size_t & count, const size_t & zpad)
  {
    // Replacement string
    std::string rs = std::to_string(count);
    if (zpad > 0) {
      std::stringstream ss;
      ss << std::setw(zpad) << std::setfill('0') << rs;
      rs = ss.str();
    }

    // Replace pattern recursively in the configuration
    util::seekAndReplace(config, pattern, rs);
  }

}  // namespace util
