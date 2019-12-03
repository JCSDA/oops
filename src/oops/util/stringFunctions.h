/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_STRINGFUNCTIONS_H_
#define OOPS_UTIL_STRINGFUNCTIONS_H_

#include <string>
#include "eckit/config/Configuration.h"

namespace util {
  namespace stringfunctions {

    void swapNameMember(const eckit::Configuration &, std::string &);

  }  // namespace stringfunctions
}  // namespace util

#endif  // OOPS_UTIL_STRINGFUNCTIONS_H_
