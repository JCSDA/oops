/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "util/abor1_cpp.h"

#include <cstdlib>
#include <string>

#include "util/Logger.h"

using oops::Log;

namespace util {
void abor1_cpp(const std::string & cderror) {
  Log::error() << cderror << std::endl;
  exit(EXIT_FAILURE);
}

void abor1_cpp(const std::string & cderror, const std::string & file,
               int line) {
  Log::error() << "ABORT: " << cderror << std::endl;
  Log::error() << "       in file '" << file << "', line " << line << std::endl;
  exit(EXIT_FAILURE);
}

void abor1_cpp_(const char cderror[]) {
  std::string s(cderror);
  util::abor1_cpp(s);
}


}  // namespace util
