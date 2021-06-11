/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_GENERALIZEDDEPARTURES_H_
#define OOPS_BASE_GENERALIZEDDEPARTURES_H_

#include "oops/util/Printable.h"

namespace oops {

/// Abstract base class for quantities
/*!
 *  The purpose of this class is to give a generic type to quantities
 *  measured by the cost function to gather them in a DualVector object.
 */

class GeneralizedDepartures : public util::Printable {
 public:
  GeneralizedDepartures() {}
  virtual ~GeneralizedDepartures() {}
 private:
  void print(std::ostream &) const = 0;
};

}  // namespace oops

#endif  // OOPS_BASE_GENERALIZEDDEPARTURES_H_
