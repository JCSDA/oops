/*
 * (C) Copyright 2017-2018  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <math.h>
#include <cfloat>

#include "oops/generic/gc99.h"

namespace oops {

double gc99(const double & distnorm) {
  // computes Gaspari-Cohn 99 localization
  // NOTE: this implementation goes to 0 at distnorm==1
  //       returns eps_double outside of distnorm>1
  // distnorm - normalized distance

  double gc99value = 0.0;

  if (distnorm < 0.5) {
    gc99value = -8.0*pow(distnorm, 5.0)+8.0*pow(distnorm, 4.0)+5.0*pow(distnorm, 3.0)-
                20.0/3.0*pow(distnorm, 2.0)+1.0;
  } else if (distnorm < 1.0) {
    gc99value = 8.0/3.0*pow(distnorm, 5.0)-8.0*pow(distnorm, 4.0)+5.0*pow(distnorm, 3.0)+
                20.0/3.0*pow(distnorm, 2.0)-10.0*distnorm+4.0-1.0/(3.0*distnorm);
  } else {
    gc99value = DBL_EPSILON;
  }
  return gc99value;
}
}  // namespace oops

