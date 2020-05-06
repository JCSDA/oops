/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <utility>

#include "eckit/config/Configuration.h"
#include "eckit/value/Value.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {
void Parameters::registerChild(ParameterBase &parameter) {
  children_.push_back(&parameter);
}
void Parameters::registerExternalParameter(std::string externalParam) {
    externalParameters_.insert(std::move(externalParam));
}
void Parameters::deserialize(const eckit::Configuration &config) {
    std::set<std::string> usedKeys = externalParameters_;
    for (auto child : children_) {
        child->deserialize(config, usedKeys);
    }
    if (spellcheckEnabled_) {
        const eckit::Value& v = config.get();
        if (!v.isNil()) {
            for (const std::string& availableKey : config.keys()) {
                if (usedKeys.find(availableKey) == usedKeys.end())
                    eckit::Log::warning() << "Warning: " << availableKey <<
                                             " was in config file, but never used." << std::endl;
            }
        }
    }
}
}  // namespace oops
