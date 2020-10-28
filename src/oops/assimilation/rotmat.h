/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_ROTMAT_H_
#define OOPS_ASSIMILATION_ROTMAT_H_

#include <cmath>

namespace oops {

/*! \file rotmat.h
 * \brief Compute the Givens rotation matrix parameters 
*/

void rotmat(const double & a, const double & b,
            double & c, double & s) {
  if (b == 0.0) {
    c = 1.0;
    s = 0.0;
  } else if (a == 0.0) {
    c = 0.0;
    s = 1.0;
  } else if (std::abs(b) > std::abs(a)) {
      double temp = a/b;
      s = 1.0/sqrt(1.0 + temp*temp);
      c = temp*s;
  } else {
      double temp = b/a;
      c = 1.0/sqrt(1.0 + temp*temp);
      s = temp*c;
  }
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_ROTMAT_H_
