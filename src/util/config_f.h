/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef UTIL_CONFIG_F_H_
#define UTIL_CONFIG_F_H_

#include<istream>

namespace eckit {
  class Configuration;
}

// -----------------------------------------------------------------------------
// These functions provide a Fortran-callable interface to Config.
// -----------------------------------------------------------------------------

namespace util {

extern "C" {
  const eckit::Configuration * config_get_testdata_f();
  bool config_element_exists_f(const eckit::Configuration & d, const char str[]);
  int config_get_data_as_int_f(const eckit::Configuration & d, const char str[]);
  double config_get_data_as_double_f(const eckit::Configuration & d, const char str[]);
  int config_get_data_length_f(const eckit::Configuration & d, const char str[]);
  void config_get_data_f(const eckit::Configuration & d, const char str[], char output[]);
}

}  // namespace util

#endif  // UTIL_CONFIG_F_H_
