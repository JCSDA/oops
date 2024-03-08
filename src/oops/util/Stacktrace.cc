/*
 * (C) Copyright 2024 The Tomorrow Companies, Inc.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef BOOST_STACKTRACE_USE_NOOP
#  include <boost/stacktrace.hpp>
#endif

#include <iostream>
#include <sstream>
#include <string>

#include "oops/util/Stacktrace.h"

namespace util {

std::string stacktrace_current() {
#ifndef BOOST_STACKTRACE_USE_NOOP
  std::ostringstream s;
  s << boost::stacktrace::stacktrace() << std::endl;
  std::string sout = s.str();
  return sout;
#else
  return "Stacktrace not available";
#endif
}

void unwind_exception_stack(const std::exception& e, std::ostream& out, int level) {
  out << "Exception: level: " << level << "\n" << e.what() << std::endl;
  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception& f) {
    unwind_exception_stack(f, out, level + 1);
  } catch (...) {
    out << "exception: level: " << level
        << "\n\tException at this level is not derived from std::exception." << std::endl;
  }
}

}  // end namespace util
