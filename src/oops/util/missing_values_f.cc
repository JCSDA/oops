/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/util/missing_values_f.h"

#include "oops/util/missingValues.h"

namespace util {

// -----------------------------------------------------------------------------
float missing_value_flt_f() {
  const float miss = util::missingValue(miss);
  return miss;
}
// -----------------------------------------------------------------------------
double missing_value_dbl_f() {
  const double miss = util::missingValue(miss);
  return miss;
}
// -----------------------------------------------------------------------------

}  // namespace util
