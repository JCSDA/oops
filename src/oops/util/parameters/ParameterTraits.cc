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

// Specialization for IntegerParameterTraits
template std::string IntegerParameterTraits<int>::valueAsJson(const int &value);

template std::string IntegerParameterTraits<size_t>::valueAsJson(const size_t &value);

template std::string IntegerParameterTraits<int64_t>::valueAsJson(const int64_t &value);

template <typename T>
std::string IntegerParameterTraits<T>::valueAsJson(const T &value) {
  return std::to_string(value);
}

// Specialization for FloatingPointParameterTraits
template std::string FloatingPointParameterTraits<float>::valueAsJson(const float &value);

template std::string FloatingPointParameterTraits<double>::valueAsJson(const double &value);

template <typename T>
std::string FloatingPointParameterTraits<T>::valueAsJson(const T &value) {
  std::stringstream valueStream;
  valueStream
    << std::setprecision(std::numeric_limits<T>::max_digits10)
    << std::scientific
    << value;
  return valueStream.str();
}

// Specialization for ParameterTraits<bool, std::false_type>
std::string ParameterTraits<bool, std::false_type>::valueAsJson(const bool &value) {
  return value ? "true" : "false";
}

// Specialization for ParameterTraits<std::string, std::false_type>
std::string ParameterTraits<std::string, std::false_type>::valueAsJson(const std::string &value) {
  return "\"" + value + "\"";
}

// Specialization for ParameterTraits<util::DateTime>
std::string ParameterTraits<util::DateTime>::valueAsJson(const util::DateTime &value) {
  return "\"" + value.toString() + "\"";
}

// Specialization for ParameterTraits<util::Duration>
std::string ParameterTraits<util::Duration>::valueAsJson(const util::Duration &value) {
  return "\"" + value.toString() + "\"";
}

// Specialization for ParameterTraits<util::PartialDateTime>
std::string ParameterTraits<util::PartialDateTime>::valueAsJson(
    const util::PartialDateTime &value)
{
  return "\"" + value.toString() + "\"";
}

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

std::string ParameterTraits<std::set<int>, std::false_type>::valueAsJson(
    const std::set<int> &value)
{
  if (value.empty()) {
    return "[]";
  }
  return "["
    + util::stringfunctions::join(
      ", ", value.begin(), value.end(), [](int n) { return std::to_string(n); })
    + "]";
}

}  // namespace oops
