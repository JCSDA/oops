/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/util/abor1_cpp.h"

#include <cstdlib>
#include <string>

#include "oops/mpi/mpi.h"
#include "oops/util/LibOOPS.h"
#include "oops/util/Logger.h"

namespace util {
void abor1_cpp(const std::string & cderror) {
  oops::Log::error() << cderror << std::endl;
  oops::LibOOPS::instance().finalise(false /* finaliseMPI? */);
  oops::mpi::world().abort(EXIT_FAILURE);
}

void abor1_cpp(const std::string & cderror, const std::string & file,
               int line) {
  oops::Log::error() << "ABORT: " << cderror << std::endl;
  oops::Log::error() << "       in file '" << file << "', line " << line << std::endl;
  oops::LibOOPS::instance().finalise(false /* finaliseMPI? */);
  oops::mpi::world().abort(EXIT_FAILURE);
}

void abor1_cpp_(const char cderror[]) {
  std::string s(cderror);
  util::abor1_cpp(s);
}


}  // namespace util
