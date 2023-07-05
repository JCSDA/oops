/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/Variables.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "eckit/types/Types.h"
#include "eckit/utils/Hash.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace oops {
// -----------------------------------------------------------------------------

Variables::Variables()
  : convention_(""), vars_(0), varMetaData_() {
  Log::trace() << "Variables::Variables" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & conf, const std::string & name)
  : convention_(""), vars_(0), channels_(0), varMetaData_() {
  Log::trace() << "Variables::Variables start " << conf << std::endl;
  std::vector<std::string> vars;
  if (!conf.get(name, vars)) {
    Log::error() << name << " not found in " << conf << std::endl;
    throw eckit::BadParameter("Undefined variable: '" + name + "'");
  }
  // hack to read channels
  if (conf.has("channels")) {
    std::string chlist = conf.getString("channels");
    std::set<int> channels = parseIntSet(chlist);
    std::copy(channels.begin(), channels.end(), std::back_inserter(channels_));
    // assuming the same channel subsetting applies to all variables
    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      for (size_t jch = 0; jch < channels_.size(); ++jch) {
        vars_.push_back(vars[jvar]+"_"+std::to_string(channels_[jch]));
      }
    }
  } else {
    vars_ = vars;
  }
  Log::trace() << "Variables::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------
Variables::Variables(const std::vector<std::string> & vars, const std::string & conv)
  : convention_(conv), vars_(vars), varMetaData_() {
  Log::trace() << "Variables::Variables start " << vars << std::endl;
  Log::trace() << "Variables::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const std::vector<std::string> & vars, const std::vector<int> & channels)
  : convention_(""), vars_(0), channels_(channels), varMetaData_() {
  Log::trace() << "Variables::Variables start " << vars << " @ " << channels << std::endl;
  if (channels.empty()) {
    vars_ = vars;
  } else {
    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      for (size_t jch = 0; jch < channels_.size(); ++jch) {
        vars_.push_back(vars[jvar]+"_"+std::to_string(channels_[jch]));
      }
    }
  }
  Log::trace() << "Variables::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & conf, const std::vector<std::string> & vars)
  : convention_(""), vars_(vars), varMetaData_(conf)
{
  Log::trace() << "Variables::Variables start " << vars << " @ " << conf << std::endl;
  Log::trace() << "Variables::Variables done" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const Variables & other)
  : convention_(other.convention_), vars_(other.vars_), channels_(other.channels_),
    varMetaData_(other.varMetaData_)
{}

// -----------------------------------------------------------------------------

Variables & Variables::operator+=(const Variables & rhs) {
  ASSERT(convention_ == rhs.convention_);
  vars_.insert(vars_.end(), rhs.vars_.begin(), rhs.vars_.end());
  // revisit late, should we add channels this way ?
  channels_.insert(channels_.end(), rhs.channels_.begin(), rhs.channels_.end());
  // remove duplicated variables and channels
  std::unordered_set<std::string> svars;
  auto mvar = std::stable_partition(vars_.begin(), vars_.end(),
        [&svars](std::string const &var) {return svars.insert(var).second;});
  vars_.erase(mvar, vars_.end());
  std::unordered_set<int> schannels;
  auto mchannel = std::stable_partition(channels_.begin(), channels_.end(),
        [&schannels](int const &channel) {return schannels.insert(channel).second;});
  channels_.erase(mchannel, channels_.end());

  // this operation adds and updates the metadata in the object from the
  // rhs object.
  for (const std::string & var : rhs.varMetaData_.keys()) {
    for (const std::string & key : rhs.varMetaData_.getSubConfiguration(var).keys()) {
      int value(-1);
      getVariableSubKeyValue(var, key, rhs.varMetaData_, value);
      setVariableSubKeyValue(var, key, value, varMetaData_);
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------

Variables & Variables::operator-=(const Variables & rhs) {
  ASSERT(convention_ == rhs.convention_);
  if (!rhs.channels().empty()) {
    throw eckit::NotImplemented(
        "Variables::operator-= not implemented for rhs objects with channels", Here());
  }
  for (auto & var : rhs.vars_) {
    vars_.erase(std::remove(vars_.begin(), vars_.end(), var), vars_.end());
  }
  return *this;
}

// -----------------------------------------------------------------------------

Variables & Variables::operator-=(const std::string & var) {
  vars_.erase(std::remove(vars_.begin(), vars_.end(), var), vars_.end());
  return *this;
}

// -----------------------------------------------------------------------------

bool Variables::operator==(const Variables & rhs) const {
  if ((convention_  != rhs.convention_) ||
      (channels_    != rhs.channels_)   ||
      (vars_.size() != rhs.vars_.size())) {
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

bool Variables::operator!=(const Variables & rhs) const {
  return (!(*this == rhs));
}

// -----------------------------------------------------------------------------

bool Variables::operator<=(const Variables & rhs) const {
  ASSERT(convention_ == rhs.convention_);
  ASSERT(channels_.empty());
  bool is_in_rhs = true;
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    is_in_rhs = is_in_rhs && rhs.has(vars_[jj]);
  }
  return is_in_rhs;
}

// -----------------------------------------------------------------------------
void Variables::addMetaData(const std::string & varname,
                            const std::string & keyname,
                            const int & keyvalue) {
  setVariableSubKeyValue(varname, keyname, keyvalue, varMetaData_);
}

// -----------------------------------------------------------------------------

void Variables::intersection(const Variables & rhs) {
  ASSERT(convention_ == rhs.convention_);
  if (this->size() * rhs.size() > 0) ASSERT(channels_.empty() == rhs.channels_.empty());

  std::vector<std::string> myvars = this->asCanonical();
  std::vector<std::string> othervars = rhs.asCanonical();
  std::vector<std::string> commonvars;
  std::set_intersection(myvars.cbegin(), myvars.cend(),
                        othervars.cbegin(), othervars.cend(), std::back_inserter(commonvars));
  vars_ = commonvars;

  if (!commonvars.empty()) {
    std::vector<int> mychannels = channels_;
    std::vector<int> otherchannels = rhs.channels_;
    std::sort(mychannels.begin(), mychannels.end());
    std::sort(otherchannels.begin(), otherchannels.end());
    std::vector<int> commonchannels;
    std::set_intersection(mychannels.cbegin(), mychannels.cend(),
                          otherchannels.cbegin(), otherchannels.cend(),
                          std::back_inserter(commonchannels));
    channels_ = commonchannels;
  } else {
    channels_.clear();
  }
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

size_t Variables::find(const std::string & var) const {
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

void Variables::push_back(const std::string & vname) {
  vars_.push_back(vname);
}

// -----------------------------------------------------------------------------

void Variables::sort() {
  std::sort(vars_.begin(), vars_.end());
  std::sort(channels_.begin(), channels_.end());
}

// -----------------------------------------------------------------------------

std::vector<std::string> Variables::asCanonical() const {
  std::vector<std::string> vars(vars_);
  std::sort(vars.begin(), vars.end());
  return vars;
}

// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  os << vars_.size() << " variables: ";
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (jj > 0) os << ", ";
    os << vars_[jj];
  }
  if (!convention_.empty()) os << " (" << convention_ << ")";

  if (!varMetaData_.empty()) os << " (" << varMetaData_ << ")";
}

// -----------------------------------------------------------------------------

int Variables::getLevels(const std::string & fieldname) const {
  int levels(-1);
  getVariableSubKeyValue(fieldname, "levels",
                         varMetaData_, levels);
  return levels;
}

// -----------------------------------------------------------------------------

void Variables::getVariableSubKeyValue(const std::string & varname,
                                       const std::string & keyname,
                                       const eckit::Configuration & variablesconf,
                                       int & intvalue) const {
  ASSERT(!variablesconf.empty());
  ASSERT(variablesconf.has(varname));
  ASSERT(variablesconf.getSubConfiguration(varname).has(keyname));
  variablesconf.getSubConfiguration(varname).get(keyname, intvalue);
}

// -----------------------------------------------------------------------------

void Variables::setVariableSubKeyValue(const std::string & varname,
                                       const std::string & keyname,
                                       const int & keyvalue,
                                       eckit::LocalConfiguration & variableslconf) {
  eckit::LocalConfiguration variablelconf =
    variableslconf.has(varname) ? variableslconf.getSubConfiguration(varname) :
                                  eckit::LocalConfiguration();
  variablelconf.set(keyname, keyvalue);
  variableslconf.set(varname, variablelconf);
}
}  // namespace oops
