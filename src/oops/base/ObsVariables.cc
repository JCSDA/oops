/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/ObsVariables.h"

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace oops {

// -----------------------------------------------------------------------------

ObsVariables::ObsVariables(const eckit::Configuration & conf, const std::string & name)
  : VariablesBase() {
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
}

// -----------------------------------------------------------------------------
ObsVariables::ObsVariables(const std::vector<std::string> & vars)
  : VariablesBase(vars), channels_({}) {
}

// -----------------------------------------------------------------------------

ObsVariables::ObsVariables(const std::vector<std::string> & vars, const std::vector<int> & channels)
  : VariablesBase(), channels_(channels) {
  if (channels.empty()) {
    vars_ = vars;
  } else {
    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      for (size_t jch = 0; jch < channels_.size(); ++jch) {
        vars_.push_back(vars[jvar]+"_"+std::to_string(channels_[jch]));
      }
    }
  }
}

// -----------------------------------------------------------------------------

ObsVariables & ObsVariables::operator+=(const ObsVariables & rhs) {
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
  return *this;
}

// -----------------------------------------------------------------------------

ObsVariables & ObsVariables::operator-=(const std::string & var) {
  vars_.erase(std::remove(vars_.begin(), vars_.end(), var), vars_.end());
  return *this;
}

// -----------------------------------------------------------------------------

bool ObsVariables::operator==(const ObsVariables & rhs) const {
  if ((channels_    != rhs.channels_)   ||
      (vars_.size() != rhs.vars_.size())) {
    return false;
  } else {
    // equivalence for the metadata is only held for the vars_list.
    std::vector<std::string> myvars = this->asCanonical();
    std::vector<std::string> othervars = rhs.asCanonical();
    return myvars == othervars;
  }
}

// -----------------------------------------------------------------------------

bool ObsVariables::operator!=(const ObsVariables & rhs) const {
  return (!(*this == rhs));
}

// -----------------------------------------------------------------------------

void ObsVariables::intersection(const ObsVariables & rhs) {
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

void ObsVariables::sort() {
  std::sort(vars_.begin(), vars_.end());
  std::sort(channels_.begin(), channels_.end());
}

// -----------------------------------------------------------------------------

}  // namespace oops
