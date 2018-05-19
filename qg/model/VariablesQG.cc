/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/VariablesQG.h"

#include<vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace qg {

// -----------------------------------------------------------------------------

VariablesQG::VariablesQG(const oops::Variables & oopsvars) {
  oops::Log::debug() << "VariablesQG oopsvar:" << oopsvars.variables() << std::endl;
  this->setF90(oopsvars.variables());
  print(oops::Log::debug());
}

// -----------------------------------------------------------------------------

VariablesQG::VariablesQG(const eckit::Configuration & config) {
  oops::Log::debug() << "VariablesQG config:" << config << std::endl;
  std::vector<std::string> vars;
  config.get("variables", vars);
  this->setF90(vars);
  print(oops::Log::debug());
}

// -----------------------------------------------------------------------------

void VariablesQG::setF90(const std::vector<std::string> vars) {
  size_t nv = vars.size();
  oops::Log::debug() << "setF90 " << nv << " vars = " << vars << std::endl;
  fvars_.resize(nv + 2);
  fvars_[0] = nv;
  for (size_t jj = 0; jj < nv; ++jj) {
     int ii = 0;
     if (vars[jj] == "x") ii = 1;
     if (vars[jj] == "q") ii = 2;
     if (vars[jj] == "u") ii = 3;
     if (vars[jj] == "v") ii = 4;
     if (vars[jj] == "bc") ii = 5;
     ASSERT(ii > 0);
     fvars_[jj+1] = ii;
  }
  fvars_[nv+1] = 999;  // just for checking
  oops::Log::debug() << "setF90 " << nv << " fvars = " << fvars_ << std::endl;
}

// -----------------------------------------------------------------------------

VariablesQG::~VariablesQG() {}

// -----------------------------------------------------------------------------

VariablesQG::VariablesQG(const VariablesQG & other): fvars_(other.fvars_) {}

// -----------------------------------------------------------------------------

void VariablesQG::print(std::ostream & os) const {
  os << "qg::VariablesQG: vars = " << fvars_;
}

// -----------------------------------------------------------------------------

}  // namespace qg
