/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_ABOR1_CPP_H_
#define OOPS_UTIL_ABOR1_CPP_H_

#include <string>

#define ABORT(s) \
  util::abor1_cpp(s, __FILE__, __LINE__)

namespace util {

/// Print an error message and abort execution
/*!
 * If called with a single string argument, the string is printed and
 * execution aborts.
 *
 * If called with 3 arguments, the first is the error message, the
 * second should be __FILE__ and the third should be __LINE__.
 * You should not call abor1_cpp directly. Use the ASSERT or ABORT macro.
 */

void abor1_cpp(const std::string & cderror);
void abor1_cpp(const std::string & cderror, const std::string & file, int line);

extern "C" {
  void abor1_cpp_(const char[]);
}

}  // namespace util

#endif  // OOPS_UTIL_ABOR1_CPP_H_
