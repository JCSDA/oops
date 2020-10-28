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

#include "oops/util/formats.h"
#include "oops/util/Printable.h"

namespace util {

// -----------------------------------------------------------------------------

class PrintAdjTest : public Printable {
 public:
  PrintAdjTest(double& dp1, double& dp2, std::string op);
  virtual ~PrintAdjTest() {}

 private:
  void print(std::ostream &) const;

  double& dp1_;
  double& dp2_;
  std::string op_;
};

// -----------------------------------------------------------------------------

PrintAdjTest::PrintAdjTest(double& dp1, double& dp2, std::string op)
  : dp1_(dp1), dp2_(dp2), op_(op)
{}

void PrintAdjTest::print(std::ostream & os) const {
  double dpdiff = std::abs(dp1_-dp2_);
  os << "  " << op_ << " adjoint test: <" << op_ << " dx1,  dx2>         : "
             << full_precision(dp1_)        << std::endl
     << "  " << op_ << " adjoint test: <dx1, " << op_ << "t dx2>         : "
             << full_precision(dp2_)        << std::endl
     << "  " << op_ << " adjoint test: difference            : "
             << full_precision(dpdiff)      << std::endl
     << "  " << op_ << " adjoint test: normalized difference : "
             << full_precision(dpdiff/dp1_) << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_PRINTADJTEST_H_
