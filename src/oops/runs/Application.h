/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_APPLICATION_H_
#define OOPS_RUNS_APPLICATION_H_

#include <iostream>
#include <string>

#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

class Application : public util::Printable {
 public:
  Application() {}
  virtual ~Application() {}
  virtual int execute(const eckit::Configuration &) const =0;

 private:
  virtual std::string appname() const =0;
  virtual void print(std::ostream & os) const {os << appname();}
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_RUNS_APPLICATION_H_
