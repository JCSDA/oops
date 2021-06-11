/*
 * (C) Copyright 2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <math.h>
#include <cfloat>

#include "oops/generic/soar.h"

namespace oops {

double soar(const double & distnorm) {
  // distnorm - normalized distance c*r, where c is the decay prameter and r is distance

  double soarvalue = (1 + distnorm)*exp(-distnorm);
  return soarvalue;
}
}  // namespace oops

