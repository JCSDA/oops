/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <iostream>
#include <iterator>
#include <utility>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/value/Value.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

void Parameters::registerChild(ParameterBase &parameter) {
  children_.push_back(&parameter);
}

void Parameters::deserialize(const eckit::Configuration &config) {
  util::CompositePath path;
  deserialize(path, config);
}

void Parameters::deserialize(util::CompositePath &path, const eckit::Configuration &config) {
  for (ParameterBase* child : children_)
    child->deserialize(path, config);
}

}  // namespace oops
