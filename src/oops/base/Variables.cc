/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/Variables.h"

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "eckit/types/Types.h"
#include "eckit/utils/Hash.h"

#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace oops {

// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & conf, const std::string & name)
  : VariablesBase(conf, name) {
}

// -----------------------------------------------------------------------------
Variables::Variables(const std::vector<std::string> & vars)
  : VariablesBase(vars) {
}

// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & conf, const std::vector<std::string> & vars)
  : VariablesBase(vars), varMetaData_(conf)
{
}

// -----------------------------------------------------------------------------

Variables & Variables::operator+=(const Variables & rhs) {
  vars_.insert(vars_.end(), rhs.vars_.begin(), rhs.vars_.end());
  // remove duplicated variables
  std::unordered_set<std::string> svars;
  auto mvar = std::stable_partition(vars_.begin(), vars_.end(),
        [&svars](std::string const &var) {return svars.insert(var).second;});
  vars_.erase(mvar, vars_.end());

  // this operation adds and updates the metadata in the object from the
  // rhs object.
  for (const std::string & var : rhs.varMetaData_.keys()) {
    eckit::LocalConfiguration varConf = rhs.varMetaData_.getSubConfiguration(var);
    varMetaData_.set(var, varConf);
  }
  return *this;
}

// -----------------------------------------------------------------------------

Variables & Variables::operator-=(const Variables & rhs) {
  for (auto & var : rhs.vars_) {
    vars_.erase(std::remove(vars_.begin(), vars_.end(), var), vars_.end());
  }
  return *this;
}

// -----------------------------------------------------------------------------

bool Variables::operator==(const Variables & rhs) const {
  if (vars_.size() != rhs.vars_.size()) {
    return false;
  } else {
    // equivalence for the metadata is only held for the vars_list.
    std::vector<std::string> myvars = this->asCanonical();
    std::vector<std::string> othervars = rhs.asCanonical();
    if (varMetaData_.empty()) {
      return myvars == othervars;
    } else {
      for (const std::string & var : vars_) {
        if ((!rhs.varMetaData_.has(var) && varMetaData_.has(var)) ||
            (rhs.varMetaData_.has(var) && !varMetaData_.has(var))) {
          return false;
        } else if (rhs.varMetaData_.has(var) && varMetaData_.has(var)) {
          std::unique_ptr<eckit::Hash> myhash(eckit::HashFactory::instance().build("md5"));
          std::unique_ptr<eckit::Hash> rhshash(eckit::HashFactory::instance().build("md5"));
          varMetaData_.getSubConfiguration(var).hash(*myhash);
          rhs.varMetaData_.getSubConfiguration(var).hash(*rhshash);
          if (!(myhash->digest() == rhshash->digest())) {
            return false;
          }
        }
      }
      return true;
    }
  }
}

// -----------------------------------------------------------------------------

Variables & Variables::operator-=(const std::string & var) {
  vars_.erase(std::remove(vars_.begin(), vars_.end(), var), vars_.end());
  return *this;
}

// -----------------------------------------------------------------------------

bool Variables::operator!=(const Variables & rhs) const {
  return (!(*this == rhs));
}

// -----------------------------------------------------------------------------

bool Variables::operator<=(const Variables & rhs) const {
  bool is_in_rhs = true;
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    is_in_rhs = is_in_rhs && rhs.has(vars_[jj]);
  }
  return is_in_rhs;
}

// -----------------------------------------------------------------------------

void Variables::intersection(const Variables & rhs) {
  std::vector<std::string> myvars = this->asCanonical();
  std::vector<std::string> othervars = rhs.asCanonical();
  std::vector<std::string> commonvars;
  std::set_intersection(myvars.cbegin(), myvars.cend(),
                        othervars.cbegin(), othervars.cend(), std::back_inserter(commonvars));
  vars_ = commonvars;
}

// -----------------------------------------------------------------------------

void Variables::sort() {
  std::sort(vars_.begin(), vars_.end());
}

// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  os << vars_.size() << " variables: ";
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (jj > 0) os << ", ";
    os << vars_[jj];
  }
  if (!varMetaData_.empty()) os << " (" << varMetaData_ << ")";
}

// -----------------------------------------------------------------------------

bool Variables::hasMetaData(const std::string & varname,
                            const std::string & keyname) const {
  bool has = false;
  if (!varMetaData_.empty()) {
    if (varMetaData_.has(varname)) {
      has = varMetaData_.getSubConfiguration(varname).has(keyname);
    }
  }
  return has;
}


// -----------------------------------------------------------------------------

int Variables::getLevels(const std::string & fieldname) const {
  int levels(-1);
  getVariableSubKeyValue(fieldname, "levels",
                         varMetaData_, levels);
  return levels;
}

}  // namespace oops
