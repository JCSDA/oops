/*
 * (C) Copyright 2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PRINTSTACKTRACE_H_
#define OOPS_UTIL_PRINTSTACKTRACE_H_

/// Print backtrace for debugging

// #define _GNU_SOURCE
// #include <boost/stacktrace.hpp>

#include "oops/util/Logger.h"

namespace oops {

void printStackTrace() {
  // Log::debug() << boost::stacktrace::stacktrace();
}

}

#endif  // OOPS_UTIL_PRINTSTACKTRACE_H_
