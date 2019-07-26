/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/VariablesQG.h"

#include<vector>

#include "eckit/config/Configuration.h"
#include "eckit/types/Types.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace qg {

// -----------------------------------------------------------------------------

VariablesQG::VariablesQG(const oops::Variables & oopsvars)
                         : oopsvars_(oopsvars) {
  oops::Log::debug() << "VariablesQG oopsvar:" << oopsvars.variables() << std::endl;
  this->setF90(oopsvars.variables());
  oops::Log::debug() << *this << std::endl;
}

// -----------------------------------------------------------------------------

VariablesQG::VariablesQG(const eckit::Configuration & config)
                         : oopsvars_(config) {
  oops::Log::debug() << "VariablesQG config:" << config << std::endl;
  std::vector<std::string> vars;
  config.get("variables", vars);
  this->setF90(vars);
  oops::Log::debug() << *this << std::endl;
}

// -----------------------------------------------------------------------------

void VariablesQG::setF90(const std::vector<std::string> vars) {
  size_t nv = vars.size();
  oops::Log::debug() << "setF90 " << nv << " vars = " << vars << std::endl;
  fvars_.resize(4);
  std::fill(fvars_.begin(), fvars_.end(), 0);
  for (size_t jj = 0; jj < nv; ++jj) {
    if (vars[jj] == "x") {
      fvars_[0] = 1;
    } else {
      if (vars[jj] == "q") {
        fvars_[1] = 1;
      } else {
        if (vars[jj] == "u") {
          fvars_[2] = 1;
        } else {
          if (vars[jj] == "v") {
            fvars_[3] = 1;
          } else {
            ABORT("Wrong variable name");
          }
        }
      }
    }
  }
  oops::Log::debug() << "setF90 " << nv << " fvars = " << fvars_ << std::endl;
}

// -----------------------------------------------------------------------------

VariablesQG::~VariablesQG() {}

// -----------------------------------------------------------------------------

VariablesQG::VariablesQG(const VariablesQG & other)
           : fvars_(other.fvars_),
             oopsvars_(other.toOopsVariables()) {}

// -----------------------------------------------------------------------------

void VariablesQG::print(std::ostream & os) const {
  os << "qg::VariablesQG: vars = " << fvars_;
}

// -----------------------------------------------------------------------------

}  // namespace qg
