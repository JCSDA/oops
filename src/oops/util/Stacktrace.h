/*
 * (C) Copyright 2024 The Tomorrow Companies, Inc.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_STACKTRACE_H_
#define OOPS_UTIL_STACKTRACE_H_

#include <exception>
#include <iostream>
#include <memory>

namespace util {

/// @brief This is a basic backport of C++23's stacktrace feature using boost.
///   The boost code is wrapped in an opaque structure that is hidden from the rest of the system.
std::string stacktrace_current();

}  // namespace util

#endif  // OOPS_UTIL_STACKTRACE_H_
