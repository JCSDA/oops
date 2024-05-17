/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/ParameterTraitsObsVariables.h"

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/ObsVariables.h"

namespace oops {

std::vector<std::string> getVariableNamesWithoutChannelSuffix(const ObsVariables &variables) {
  if (variables.channels().empty())
    return variables.variables();

  auto throwException = [&variables] {
    std::stringstream msg;
    msg << "The ObsVariables object " << variables << " cannot be stored in a Configuration "
           "object: not all variables have the same channels";
    throw std::runtime_error(msg.str());
  };

  if (variables.variables().size() % variables.channels().size() != 0)
    throwException();  // clearly, not all variables have the same channel suffixes

  typedef std::set<std::string> ChannelSuffixes;
  ChannelSuffixes expectedChannels;
  for (int channel : variables.channels())
    expectedChannels.insert(std::to_string(channel));

  // Split each variable name into a base name and a channel suffix
  const char channelSeparator = '_';
  std::vector<std::string> uniqueBaseNames;
  std::map<std::string, ChannelSuffixes> channelsPerVariable;
  for (const std::string& name : variables.variables()) {
    const std::string::size_type separatorPosition = name.find_last_of(channelSeparator);
    if (separatorPosition == std::string::npos)
      throwException();  // no channel suffix
    std::string baseName = name.substr(0, separatorPosition);
    std::string channel = name.substr(separatorPosition + 1);

    std::map<std::string, ChannelSuffixes>::iterator it = channelsPerVariable.find(baseName);
    if (it == channelsPerVariable.end()) {
      // It's the first time we're encountering this base name
      uniqueBaseNames.push_back(baseName);
      it = channelsPerVariable.insert(std::make_pair(baseName, ChannelSuffixes())).first;
    }
    it->second.insert(std::move(channel));
  }

  if (std::any_of(channelsPerVariable.begin(), channelsPerVariable.end(),
                  [&expectedChannels]
                  (const std::pair<const std::string, ChannelSuffixes> &nameAndChannels)
                  { return nameAndChannels.second != expectedChannels; }))
    throwException();  // some variables have different channel suffixes than expected

  return uniqueBaseNames;
}

std::string ParameterTraits<ObsVariables>::valueAsJson(const ObsVariables &value) {
  const std::vector<std::string> varNames = getVariableNamesWithoutChannelSuffix(value);
  if (varNames.empty()) {
    return "[]";
  }
  return "["
    + util::stringfunctions::join(
      ", ", varNames.begin(), varNames.end(), [](std::string s) { return "\"" + s + "\""; })
    + "]";
}

}  // namespace oops
