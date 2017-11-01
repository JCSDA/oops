/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cstring>
#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "util/config_f.h"

using oops::Log;

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
    conf_->set("grandparent.son.maxint", long(imax));
    conf_->set("grandparent.son.minint", long(imin));
    conf_->set("grandparent.son.maxdouble", dmax);
    conf_->set("grandparent.son.mindouble", dmin);
    conf_->set("grandparent.daughter.grandson", long(21));
    conf_->set("grandparent.daughter.hamster", "Errol");
  }

  ~ConfigTestData() {}

  boost::scoped_ptr<eckit::LocalConfiguration> conf_;
};

// -----------------------------------------------------------------------------

const eckit::Configuration * config_get_testdata_f() {
  const eckit::Configuration * cfg = &ConfigTestData::getConfig();
  return cfg;
}

bool config_element_exists_f(const eckit::Configuration & d, const char str[]) {
  const std::string s(str);
  Log::debug() << "config_element_exists_f #" << s << "#" << std::endl;
  return d.has(s);
}

int config_get_data_as_int_f(const eckit::Configuration & d, const char str[]) {
  const std::string s(str);
  Log::debug() << "config_get_data_as_int_f #" << s << "# from " << d << std::endl;
  std::istringstream is;
  is.str(d.getString(s));
  int i;
  is >> i;
  return i;
}

double config_get_data_as_double_f(const eckit::Configuration & d, const char str[]) {
  const std::string s(str);
  Log::debug() << "config_get_data_as_double_f #" << s << "#" << std::endl;
  std::istringstream is;
  is.str(d.getString(s));
  double dd;
  is >> dd;
  return dd;
}

int config_get_data_length_f(const eckit::Configuration & d, const char str[]) {
  const std::string s(str);
  Log::debug() << "config_get_data_length_f #" << s << "#" << std::endl;
  const std::string data(d.getString(s));
  return data.length();
}

void config_get_data_f(const eckit::Configuration & d, const char str[], char output[]) {
  const std::string s(str);
  Log::debug() << "config_get_data_f #" << s << "#" << std::endl;
  const std::string data(d.getString(s));
  size_t lendata = data.size();
  strncpy(output, data.c_str(), lendata);
}

// -----------------------------------------------------------------------------

}  // namespace util
