/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_PRINTABLE_H_
#define OOPS_UTIL_PRINTABLE_H_

#include <ostream>

namespace util {

// -----------------------------------------------------------------------------

class Printable {
 public:
  Printable() {}
  virtual ~Printable() {}

/// Print human readable informations
  friend std::ostream & operator<< (std::ostream & os, const Printable & self) {
    self.print(os);
    return os;
  }

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_PRINTABLE_H_
