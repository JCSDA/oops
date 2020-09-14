/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/parameters/ConfigurationParameter.h"

namespace oops {

void ConfigurationParameter::deserialize(util::CompositePath &path,
                                         const eckit::Configuration &config) {
  value_ = eckit::LocalConfiguration(config);
}

void ConfigurationParameter::serialize(eckit::LocalConfiguration &config) const {
  for (const std::string& key : value_.keys())
    config.set(key, value_.getSubConfiguration(key));
}

}  // namespace oops
