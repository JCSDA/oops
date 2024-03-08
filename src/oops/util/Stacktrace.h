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
#include <string>

namespace util {

/// @brief This is a basic backport of C++23's stacktrace feature.
/// @details This uses Boost instead of eckit's Backtrace code because
///   the eckit implementation has a lot of platform-specific code
///   and libstdc++-related code for symbol demangling that is more
///   generalized by Boost. Boost's stacktrace code also provides
///   line numbers whenever possible.
/// @returns A string with the current stacktrace.
std::string stacktrace_current();

/// @brief Print details of the exception stack.
/// @details This C++11 feature was added after eckit's Exception classes
///   were written. It is more general, as all of the JEDI-related exceptions
///   ultimately derive from std::exception, including the eckit-related ones.
/// @param e is the exception.
/// @param out is the output stream.
/// @param level denotes current depth in the stack. Used because
///   unwind_exception_stack is a recursive function.
void unwind_exception_stack(const std::exception& e, std::ostream& out = std::cerr, int level = 0);

}  // namespace util

#endif  // OOPS_UTIL_STACKTRACE_H_
