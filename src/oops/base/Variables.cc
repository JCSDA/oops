/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/base/Variables.h"

#include <iostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "util/Logger.h"
#include "util/abor1_cpp.h"

// -----------------------------------------------------------------------------
namespace oops {
// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & conf)
 : convention_(""), vars_(0), conf_(conf)
{
  std::string vars, s;
  conf.get("variables", vars);

  std::istringstream svars(vars);
  while (std::getline(svars, s, ',' ))
    vars_.push_back(s);
}

// -----------------------------------------------------------------------------

Variables::Variables(const std::vector<std::string> & vars, const std::string & conv)
 : convention_(conv), vars_(vars), conf_()
{
  std::string svars = "";
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (jj > 0) svars += ", ";
    svars += vars_[jj];
  }
  eckit::LocalConfiguration conf;
  conf_.set("nvars", vars_.size());
  conf_.set("variables", svars);
}

// -----------------------------------------------------------------------------

Variables::Variables(const Variables & other)
 : convention_(other.convention_), vars_(other.vars_), conf_(other.conf_)
{}

// -----------------------------------------------------------------------------

Variables::~Variables() {}

// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  os << vars_.size() << " variables: ";
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (jj > 0) os << ", ";
    os << vars_[jj];
  }
  if (!convention_.empty()) os << " (" << convention_ << ")";
}

// -----------------------------------------------------------------------------

}  // namespace oops
