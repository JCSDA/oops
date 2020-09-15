/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_INTSETPARSER_H_
#define OOPS_UTIL_INTSETPARSER_H_

#include <set>
#include <string>

// Import the definition of contains(), which used to be declared here
#include "oops/util/AssociativeContainers.h"

namespace oops {
  std::set<int> parseIntSet(const std::string &);
}  // namespace oops

#endif  // OOPS_UTIL_INTSETPARSER_H_
