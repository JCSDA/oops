/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/VariablesBase.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "eckit/types/Types.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace oops {

// -----------------------------------------------------------------------------

VariablesBase::VariablesBase(const eckit::Configuration & conf, const std::string & name) {
  if (!conf.get(name, vars_)) {
    Log::error() << name << " not found in " << conf << std::endl;
    throw eckit::BadParameter("Undefined variable: '" + name + "'");
  }
}

// -----------------------------------------------------------------------------
VariablesBase::VariablesBase(const std::vector<std::string> & vars)
  : vars_(vars) {
}

// -----------------------------------------------------------------------------

bool VariablesBase::has(const std::string & var) const {
  bool found = false;
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    found = found || vars_[jj] == var;
  }
  return found;
}

// -----------------------------------------------------------------------------

size_t VariablesBase::find(const std::string & var) const {
  size_t ii = vars_.size();
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (vars_[jj] == var) ii = jj;
  }
  if (ii >= vars_.size()) {
    Log::error() << "Could not find " << var << " in variables list: " << vars_ << std::endl;
  }
  ASSERT(ii < vars_.size());
  return ii;
}

// -----------------------------------------------------------------------------

void VariablesBase::push_back(const std::string & vname) {
  vars_.push_back(vname);
}

// -----------------------------------------------------------------------------

std::vector<std::string> VariablesBase::asCanonical() const {
  std::vector<std::string> vars(vars_);
  std::sort(vars.begin(), vars.end());
  return vars;
}

// -----------------------------------------------------------------------------

void VariablesBase::print(std::ostream & os) const {
  os << vars_.size() << " Variables: ";
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (jj > 0) os << ", ";
    os << vars_[jj];
  }
}

}  // namespace oops
