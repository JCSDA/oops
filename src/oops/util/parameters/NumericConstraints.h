/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARAMETERS_NUMERICCONSTRAINTS_H_
#define OOPS_UTIL_PARAMETERS_NUMERICCONSTRAINTS_H_

#include <exception>
#include <memory>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/ParameterConstraint.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Constrains a numeric parameter of type T to be greater than or equal to a specified
/// value.
template <typename T>
class MinConstraint : public ParameterConstraint<T> {
 public:
  explicit MinConstraint(const T &min) : min_(min) {}

  void checkValue(const std::string &path, const T &value) const override {
    if (!(value >= min_))
      throw eckit::BadValue(path + ": Property value " + std::to_string(value) +
                            " is not greater than or equal to " + std::to_string(min_));
  }

  PropertyJsonSchema jsonSchema() const override {
    return {{"minimum", std::to_string(min_)}};
  }

 private:
  T min_;
};

template <typename T>
std::shared_ptr<MinConstraint<T>> minConstraint(const T &min) {
  return std::make_shared<MinConstraint<T>>(min);
}

/// \brief Constrains a numeric parameter of type T to be greater than a specified value.
template <typename T>
class ExclusiveMinConstraint : public ParameterConstraint<T> {
 public:
  explicit ExclusiveMinConstraint(const T &min) : min_(min) {}

  void checkValue(const std::string &path, const T &value) const override {
    if (!(value > min_))
      throw eckit::BadValue(path + ": Property value " + std::to_string(value) +
                            " is not greater than " + std::to_string(min_));
  }

  PropertyJsonSchema jsonSchema() const override {
    return {{"exclusiveMinimum", std::to_string(min_)}};
  }

 private:
  T min_;
};

template <typename T>
std::shared_ptr<ExclusiveMinConstraint<T>> exclusiveMinConstraint(const T &min) {
  return std::make_shared<ExclusiveMinConstraint<T>>(min);
}

/// \brief Constrains a numeric parameter of type T to be less than or equal to a specified
/// value.
template <typename T>
class MaxConstraint : public ParameterConstraint<T> {
 public:
  explicit MaxConstraint(const T &max) : max_(max) {}

  void checkValue(const std::string &path, const T &value) const override {
    if (!(value <= max_))
      throw eckit::BadValue(path + ": Property value " + std::to_string(value) +
                            " is not less than or equal to " + std::to_string(max_));
  }

  PropertyJsonSchema jsonSchema() const override {
    return {{"maximum", std::to_string(max_)}};
  }

 private:
  T max_;
};

template <typename T>
std::shared_ptr<MaxConstraint<T>> maxConstraint(const T &max) {
  return std::make_shared<MaxConstraint<T>>(max);
}

/// \brief Constrains a numeric parameter of type T to be less than a specified value.
template <typename T>
class ExclusiveMaxConstraint : public ParameterConstraint<T> {
 public:
  explicit ExclusiveMaxConstraint(const T &max) : max_(max) {}

  void checkValue(const std::string &path, const T &value) const override {
    if (!(value < max_))
      throw eckit::BadValue(path + ": Property value " + std::to_string(value) +
                            " is not less than " + std::to_string(max_));
  }

  PropertyJsonSchema jsonSchema() const override {
    return {{"exclusiveMaximum", std::to_string(max_)}};
  }

 private:
  T max_;
};

template <typename T>
std::shared_ptr<ExclusiveMaxConstraint<T>> exclusiveMaxConstraint(const T &max) {
  return std::make_shared<ExclusiveMaxConstraint<T>>(max);
}

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_NUMERICCONSTRAINTS_H_
