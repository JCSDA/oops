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
 : convention_(""), vars_(0), fvars_(0)
{
  conf.get("variables", vars_);
  this->setF90();
}

// -----------------------------------------------------------------------------

Variables::Variables(const std::vector<std::string> & vars, const std::string & conv)
 : convention_(conv), vars_(vars), fvars_(0)
{
  this->setF90();
}

// -----------------------------------------------------------------------------

Variables::Variables(const Variables & other)
 : convention_(other.convention_), vars_(other.vars_), fvars_(other.fvars_)
{}

// -----------------------------------------------------------------------------

Variables::~Variables() {}

// -----------------------------------------------------------------------------

// This should be re-written to depend on naming convention
void Variables::setF90() {
  if (!convention_.empty()) ABORT("Variables::setF90 not implemented");
  size_t nv = vars_.size();
  Log::debug() << "setF90 " << nv << " vars = " << vars_ << std::endl;
  fvars_.resize(nv + 2);
  fvars_[0] = nv;
  for (size_t jj = 0; jj < nv; ++jj) {
     int ii = 0;
     if (vars_[jj]=="x") ii = 1;
     if (vars_[jj]=="q") ii = 2;
     if (vars_[jj]=="u") ii = 3;
     if (vars_[jj]=="v") ii = 4;
     if (vars_[jj]=="bc") ii = 5;
//     ASSERT(ii > 0);
     fvars_[jj+1] = ii;
  }
  fvars_[nv+1] = 999;  // just for checking
  Log::debug() << "setF90 " << nv << " fvars = " << fvars_ << std::endl;
}

// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (jj > 0) os << ", ";
    os << vars_[jj];
  }
  if (!convention_.empty()) os << " (" << convention_ << ")";
}

// -----------------------------------------------------------------------------

}  // namespace oops
