/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/base/Variables.h"

#include <iostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace oops {
// -----------------------------------------------------------------------------

Variables::Variables()
  : convention_(""), vars_(0), conf_(), fconf_() {
  Log::trace() << "Variables::Variables" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & conf)
  : convention_(""), vars_(0), conf_(), fconf_() {
  Log::trace() << "Variables::Variables start " << conf << std::endl;
  conf.get("variables", vars_);
  this->setConf();
  Log::trace() << "Variables::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const std::vector<std::string> & vars, const std::string & conv)
  : convention_(conv), vars_(vars), conf_(), fconf_() {
  Log::trace() << "Variables::Variables start " << vars << std::endl;
  this->setConf();
  Log::trace() << "Variables::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const Variables & other)
  : convention_(other.convention_), vars_(other.vars_), conf_(other.conf_), fconf_(other.fconf_)
{}

// -----------------------------------------------------------------------------

Variables & Variables::operator+=(const Variables & rhs) {
  ASSERT(convention_ == rhs.convention_);
  vars_.insert(vars_.end(), rhs.vars_.begin(), rhs.vars_.end());
  this->setConf();
  return *this;
}

// -----------------------------------------------------------------------------

bool Variables::has(const std::string & var) const {
  bool found = false;
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    found = found || vars_[jj] == var;
  }
  return found;
}

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

void Variables::setConf() {
  conf_.set("nvars", vars_.size());
  conf_.set("variables", vars_);

  std::string svars = "";
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (jj > 0) svars += ", ";
    svars += vars_[jj];
  }
  fconf_.set("nvars", vars_.size());
  fconf_.set("variables", svars);
}

// -----------------------------------------------------------------------------

}  // namespace oops
