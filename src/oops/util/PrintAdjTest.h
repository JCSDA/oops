/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_PRINTADJTEST_H_
#define OOPS_UTIL_PRINTADJTEST_H_

/// Wraps printing out of the result of an adjoint test.

#include <iostream>
#include <string>

#include "oops/util/Printable.h"

namespace util {

// -----------------------------------------------------------------------------

class PrintAdjTest : public Printable {
 public:
  PrintAdjTest(const double& dp1, const double& dp2, std::string op);
  virtual ~PrintAdjTest() {}

 private:
  void print(std::ostream &) const;

  const double& dp1_;
  const double& dp2_;
  std::string op_;
};

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_PRINTADJTEST_H_
