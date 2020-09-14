/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERS_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERS_H_

#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "oops/util/parameters/ParameterBase.h"

namespace eckit {
  class Channel;
  class Configuration;
}

namespace oops {

/// \brief This macro needs to be invoked at the top of the declaration of each abstract subclass
/// of Parameters, with \p className set to the name of the class being declared and \p
/// baseClassName set to the name of the class from which it (directly) inherits.
///
/// It implements constructors and assignment operators in such a way that
/// * the children_ vector is set up correctly (for example, after copy construction it contains
///   pointers to members the new Parameters instance rather than the instance that was copied);
/// * the abstract subclass provides a clone() method;
/// * the abstract subclass cannot be instantiated (the constructors are protected);
/// * the assignment operators are protected to prevent slicing.
#define OOPS_ABSTRACT_PARAMETERS(className, baseClassName) \
 protected: \
  className() : baseClassName() {} \
  explicit className(Parameters* parent) : baseClassName(parent) {} \
  className(const className &other) : className() { *this = other; } \
  className(className &&other) : className() { *this = std::move(other); } \
  className &operator=(const className &) = default; \
  className &operator=(className &&) = default; \
 private: \
  virtual className* cloneImpl() const = 0; \
 public: \
  std::unique_ptr<className> clone() const { \
    return std::unique_ptr<className>(cloneImpl()); \
  } \
 private:

/// \brief This macro needs to be invoked at the top of the declaration of each concrete subclass
/// of Parameters, with \p className set to the name of the class being declared and \p
/// baseClassName set to the name of the class from which it (directly) inherits.
///
/// It implements constructors and assignment operators in such a way that
/// * the class can be instantiated using a constructor taking an optional parameter
///   `Parameters* parent = nullptr`. If `parent` is non-null, the new instance will be registered
///   as its child, so that deserialization of `parent` will trigger deserialization of the new
///   instance;
/// * the children_ vector is set up correctly (for example, after copy construction it contains
///   pointers to members the new Parameters instance rather than the instance that was copied);
/// * the concrete subclass provides a clone() method.
#define OOPS_CONCRETE_PARAMETERS(className, baseClassName) \
 public: \
  className() : baseClassName() {} \
  explicit className(oops::Parameters* parent) : baseClassName(parent) {} \
  className(const className &other) : className() { *this = other; } \
  className(className &&other) : className() { *this = std::move(other); } \
  className &operator=(const className &) = default; \
  className &operator=(className &&) = default; \
 private: \
  virtual className* cloneImpl() const { \
    return new className(*this); \
  } \
 public: \
  std::unique_ptr<className> clone() const { \
    return std::unique_ptr<className>(cloneImpl()); \
  } \
 private:

/// \brief Abstract base class for collections of parameters located at the same level of the
/// parameter hierarchy.
///
/// Abstract subclasses should invoke the OOPS_ABSTRACT_PARAMETERS macro at the top of
/// the class definition. Concrete subclasses should invoke the OOPS_CONCRETE_PARAMETERS
/// macro instead. Example:
///
///     class MyAbstractParameters : public Parameters {
///       OOPS_ABSTRACT_PARAMETERS(MyAbstractParameters, Parameters)
///       // rest of class definition
///     };
///
///     class MyConcreteParameters : public MyAbstractParameters {
///       OOPS_CONCRETE_PARAMETERS(MyConcreteParameters, MyAbstractParameters)
///       // rest of class definition
///     };
class Parameters : public ParameterBase {
 protected:
  Parameters() : ParameterBase() {}
  explicit Parameters(Parameters* parent) : ParameterBase(parent) {}
  /// \brief Copy constructor. Doesn't copy the list of children.
  Parameters(const Parameters &) : Parameters() {}
  /// \brief Move constructor. Doesn't move the list of children.
  Parameters(Parameters &&) : Parameters() {}
  /// \brief Copy assignment operator. Doesn't assign the list of children.
  Parameters &operator=(const Parameters &) { return *this; }
  /// \brief Move assignment operator. Doesn't assign the list of children.
  Parameters &operator=(Parameters &&) { return *this; }

 private:
  virtual Parameters* cloneImpl() const = 0;

 public:
  /// \brief Return a unique_ptr to a clone of this object.
  std::unique_ptr<Parameters> clone() const {
    return std::unique_ptr<Parameters>(cloneImpl());
  }

  /// \brief Add \p parameter to the list of parameters processed by subsequent calls to
  /// deserialize().
  void registerChild(ParameterBase &parameter);

  // Import the version of deserialize declared in the base class. We will override it and
  // add an extra overload.
  using ParameterBase::deserialize;

  /// \brief Load the values of all previously registered parameters from the top-level
  /// configuration \p config.
  void deserialize(const eckit::Configuration &config);

  /// \brief Load the values of all previously registered parameters from the (not necessarily
  /// top-level) configuration \p config.
  ///
  /// \param path
  ///   Location of the configuration \p config in the full configuration loaded from a YAML file.
  ///   This object is modified in-place during deserialization, but restored to its original state
  ///   before this function returns.
  /// \param config
  ///   Configuration from which parameter values are to be loaded.
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

 private:
  std::vector<ParameterBase*> children_;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERS_H_
