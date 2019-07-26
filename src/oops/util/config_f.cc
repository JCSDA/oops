/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/util/config_f.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

namespace util {

// -----------------------------------------------------------------------------

class ConfigTestData {
 public:
  static const eckit::Configuration & getConfig() {return *getInstance().conf_;}

 private:
  static ConfigTestData& getInstance() {
    static ConfigTestData theData;
    return theData;
  }

  ConfigTestData() : conf_(new eckit::LocalConfiguration()) {
    int imin = std::numeric_limits<int>::min();
    int imax = std::numeric_limits<int>::max();
    double dmax = std::numeric_limits<double>::max();
    double dmin = -dmax;

    conf_->set("grandparent.son.granddaughter", 3.14159);
    conf_->set("grandparent.son.gerbil", "");
    conf_->set("grandparent.son.maxint", imax);
    conf_->set("grandparent.son.minint", imin);
    conf_->set("grandparent.son.maxdouble", dmax);
    conf_->set("grandparent.son.mindouble", dmin);
    conf_->set("grandparent.daughter.grandson", 21);
    conf_->set("grandparent.daughter.hamster", "Errol");
  }

  ~ConfigTestData() {}

  std::unique_ptr<eckit::LocalConfiguration> conf_;
};

// -----------------------------------------------------------------------------

const eckit::Configuration * config_get_testdata_f() {
  const eckit::Configuration * cfg = &ConfigTestData::getConfig();
  return cfg;
}

bool config_element_exists_f(const eckit::Configuration & d, const char str[]) {
  const std::string s(str);
  return d.has(s);
}

int config_get_data_as_int_f(const eckit::Configuration & d, const char str[]) {
  const std::string s(str);
  return d.getInt(s);
}

double config_get_data_as_double_f(const eckit::Configuration & d, const char str[]) {
  const std::string s(str);
  return d.getDouble(s);
}

int config_get_data_length_f(const eckit::Configuration & d, const char str[]) {
  const std::string s(str);
  const std::string data(d.getString(s));
  return data.length();
}

void config_get_data_f(const eckit::Configuration & d, const char str[], char output[]) {
  const std::string s(str);
  const std::string data(d.getString(s));
  size_t lendata = data.size();
  strncpy(output, data.c_str(), lendata);
}

int config_get_data_dimension_f(const eckit::Configuration & d, const char str[]) {
  const std::string s(str);
  const std::vector < std::string > data(d.getStringVector(s));
  return data.size();
}

int config_get_data_element_length_f(const eckit::Configuration & d, const char str[],
                                     const int & index) {
  const std::string s(str);
  const std::vector < std::string > data(d.getStringVector(s));
  return data[index-1].length();
}

void config_get_data_element_f(const eckit::Configuration & d, const char str[],
                               const int & index, char output[]) {
  const std::string s(str);
  const std::vector < std::string > data(d.getStringVector(s));
  size_t lendata = data[index-1].size();
  strncpy(output, data[index-1].c_str(), lendata);
}

/*! Vector retrival from config file
 * 
 * 
 * \warning It is the responsibility of the calling Fortran routine to allocate these
 *
 */ 

  void config_get_double_vector_f(const eckit::Configuration & d, const char str[],
                                  const std::size_t & length, double * vec) {
  const std::string s(str);
  std::vector<double> x = d.getDoubleVector(s);
  ASSERT(x.size() == length);
  for (std::size_t jj = 0; jj < length; ++jj) vec[jj] = x[jj];
}

  void config_get_float_vector_f(const eckit::Configuration & d, const char str[],
                                 const std::size_t & length, float * vec) {
  const std::string s(str);
  std::vector<float> x = d.getFloatVector(s);
  ASSERT(x.size() == length);
  for (std::size_t jj = 0; jj < length; ++jj) vec[jj] = x[jj];
}

  void config_get_int_vector_f(const eckit::Configuration & d, const char str[],
                               const std::size_t & length, std::int32_t * vec) {
  const std::string s(str);
  std::vector<int> x = d.getIntVector(s);
  ASSERT(x.size() == length);
  for (std::size_t jj = 0; jj < length; ++jj) vec[jj] = x[jj];
}

// -----------------------------------------------------------------------------

}  // namespace util
