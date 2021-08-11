/*
 * (C) Copyright 2021 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/util/parameters/OptionalParameter.h"

#include "eckit/value/Value.h"
#include "oops/util/parameters/ObjectJsonSchema.h"

namespace oops {

namespace {

/// Same as eckit::LocalConfiguration, except that it offers an extra public constructor
/// taking an eckit::Value.
///
/// This is needed to implement the OptionalParameter<void>::serialize() method; it doesn't seem
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

void OptionalParameter<void>::deserialize(util::CompositePath &/*path*/,
                                          const eckit::Configuration &config) {
  value_ = config.has(name_);
}

void OptionalParameter<void>::serialize(eckit::LocalConfiguration &config) const {
  if (value_)
    config.set(name_, LocalConfigurationEx(eckit::Value()));
}

ObjectJsonSchema OptionalParameter<void>::jsonSchema() const {
  ObjectJsonSchema schema = ObjectJsonSchema({{name_, {{"type", "\"null\""}}}});
  if (description_ != "") {
    schema.extendPropertySchema(name_, {{"description", "\"" + description_ + "\""}});
  }
  return schema;
}

}  // namespace oops
