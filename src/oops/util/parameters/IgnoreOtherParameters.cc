/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/parameters/IgnoreOtherParameters.h"
#include "oops/util/parameters/ObjectJsonSchema.h"

namespace oops {

ObjectJsonSchema IgnoreOtherParameters::jsonSchema() const {
  // An empty schema, imposing no constraints.
  return ObjectJsonSchema({}, {}, true /*additionalProperties?*/);
}

}  // namespace oops
