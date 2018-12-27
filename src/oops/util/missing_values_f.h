/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_MISSING_VALUES_F_H_
#define OOPS_UTIL_MISSING_VALUES_F_H_

namespace util {

extern "C" {
  float missing_value_flt_f();
  double missing_value_dbl_f();
}

}  // namespace util

#endif  // OOPS_UTIL_MISSING_VALUES_F_H_
