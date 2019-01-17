/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_CONFIG_F_H_
#define OOPS_UTIL_CONFIG_F_H_

#include<istream>
#include<cstdint>

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
  int config_get_data_dimension_f(const eckit::Configuration & d, const char str[]);
  int config_get_data_element_length_f(const eckit::Configuration & d, const char str[],
                                       const int & index);
  void config_get_data_element_f(const eckit::Configuration & d, const char str[],
                                 const int & index, char output[]);
  void config_get_double_vector_f(const eckit::Configuration & d, const char str[],
                                        const std::size_t & length, double * vec);
  void config_get_float_vector_f(const eckit::Configuration & d, const char str[],
                                        const std::size_t & length, float * vec);
  void config_get_int_vector_f(const eckit::Configuration & d, const char str[],
			       const std::size_t & length, std::int32_t * vec);
}

}  // namespace util

#endif  // OOPS_UTIL_CONFIG_F_H_
