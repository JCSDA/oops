
/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/ConfigHelpers.h"

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/ConfigFunctions.h"

namespace util {

void setMember(eckit::LocalConfiguration & conf, const int mem) {
  if (conf.has("member pattern")) {
    const std::string memberPattern = conf.getString("member pattern");
    util::seekAndReplace(conf, memberPattern, std::to_string(mem));
  } else {
    conf.set("member", mem);
  }
}

}  // namespace util
