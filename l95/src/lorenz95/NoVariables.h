/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_NOVARIABLES_H_
#define LORENZ95_NOVARIABLES_H_

#include <iostream>

#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {

class NoVariables : public util::Printable {
 public:
  NoVariables() {}
  explicit NoVariables(const eckit::Configuration &) {}
  ~NoVariables() {}

 private:
  void print(std::ostream & os) const {os << "NoVariables";}
};

}  // namespace lorenz95

#endif  // LORENZ95_NOVARIABLES_H_
