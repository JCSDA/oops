/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/parameters/ParameterTraits.h"

#include "eckit/types/Types.h"
#include "eckit/value/Value.h"
#include "oops/util/IntSetParser.h"

namespace oops {

namespace {

/// Same as eckit::LocalConfiguration, except that it offers an extra public constructor
/// taking an eckit::Value.
///
/// This is needed to implement the ParameterTraits<util::Null>::set() method; it doesn't seem
/// possible to set the value of a Configuration option to NilValue by calling any public method of
/// LocalConfiguration.
///
/// This is obviously a hack, hence we declare this class here rather than in a separate header
/// so as not to encourage its wider use.
class LocalConfigurationEx : public eckit::LocalConfiguration {
 public:
  // This class supports constructors taking the same parameters as the public constructors
  // of LocalConfiguration...
  explicit LocalConfigurationEx(char separator = '.') : LocalConfiguration(separator) {}
  explicit LocalConfigurationEx(eckit::Stream& s) : LocalConfiguration(s) {}
  explicit LocalConfigurationEx(const eckit::Configuration& other) : LocalConfiguration(other) {}
  LocalConfigurationEx(const Configuration& other, const std::string& path)
    : LocalConfiguration(other, path) {}

  // and in addition this constructor, which calls a protected constructor of
  // LocalConfiguration.
  explicit LocalConfigurationEx(const eckit::Value& value, char separator = '.')
    : LocalConfiguration(value, separator) {}
};

}  // namespace

// Specialization for std::set<int>
boost::optional<std::set<int>> ParameterTraits<std::set<int>, std::false_type>::get(
    util::CompositePath &path, const eckit::Configuration &config, const std::string &name) {
  boost::optional<std::set<int>> result;

  if (!config.has(name))
    return result;

  const std::string valueAsString = config.getString(name);
  try {
    result = oops::parseIntSet(valueAsString);
  } catch (eckit::Exception &) {
    throw eckit::Exception(path.path() + ": '" + valueAsString +
                           "' isn't a list of comma-separated integers or ranges of integers",
                           Here());
  }

  return result;
}

void ParameterTraits<std::set<int>, std::false_type>::set(
    eckit::LocalConfiguration &config, const std::string &name, const std::set<int> &value) {
  std::stringstream stream;
  stream << value;
  std::string valueAsString = stream.str();
  valueAsString.erase(0, 1);  // erase the first character (opening bracket)
  valueAsString.erase(valueAsString.size() - 1, 1);  // erase the last character (closing bracket)
  config.set(name, valueAsString);
}

ObjectJsonSchema ParameterTraits<std::set<int>, std::false_type>::jsonSchema(
    const std::string &name) {
  return ParameterTraits<std::string>::jsonSchema(name);
}

}  // namespace oops
