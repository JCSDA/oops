/*
 * (C) Copyright 2024 The Tomorrow Companies, Inc.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef BOOST_STACKTRACE_USE_NOOP
#  include <boost/stacktrace.hpp>
#endif

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

}  // end namespace util
