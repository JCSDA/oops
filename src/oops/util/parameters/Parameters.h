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
#include "oops/util/Printable.h"

namespace eckit {
  class Channel;
  class Configuration;
}

namespace oops {

/// \brief This macro may be invoked at the top of the declaration of an abstract subclass of
/// Parameters instead of OOPS_ABSTRACT_PARAMETERS() if the definition of the (non-copy and
/// non-move) constructors needs to be customized. See OOPS_ABSTRACT_PARAMETERS() for more
/// information.
#define OOPS_ABSTRACT_PARAMETERS_ENABLE_COPY_AND_MOVE(className, baseClassName) \
 protected: \
  className(const className &other) : className() { *this = other; } \
  className(className &&other) : className() { *this = std::move(other); } \
  className &operator=(const className &) = default; \
  className &operator=(className &&) = default; \
 private: \
  className* cloneImpl() const override = 0; \
 public: \
  std::unique_ptr<className> clone() const { \
    return std::unique_ptr<className>(cloneImpl()); \
  } \
 private:

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
///
/// If you need to customize the definition of the (non-copy and non-move) constructors, call
/// `OOPS_ABSTRACT_PARAMETERS_ENABLE_COPY_AND_MOVE(className, baseClassName)` instead of this macro
/// and define the constructors on your own.
#define OOPS_ABSTRACT_PARAMETERS(className, baseClassName) \
 protected: \
  className() : baseClassName() {} \
  explicit className(Parameters* parent) : baseClassName(parent) {} \
  OOPS_ABSTRACT_PARAMETERS_ENABLE_COPY_AND_MOVE(className, baseClassName)

/// \brief This macro may be invoked at the top of the declaration of a concrete subclass of
/// Parameters instead of OOPS_CONCRETE_PARAMETERS() if the definition of the (non-copy and
/// non-move) constructors needs to be customized. See OOPS_CONCRETE_PARAMETERS() for more
/// information.
#define OOPS_CONCRETE_PARAMETERS_ENABLE_COPY_AND_MOVE(className, baseClassName) \
 public: \
  className(const className &other) : className() { *this = other; } \
  className(className &&other) : className() { *this = std::move(other); } \
  className &operator=(const className &) = default; \
  className &operator=(className &&) = default; \
 private: \
  className* cloneImpl() const override { \
    return new className(*this); \
  } \
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
///
/// If you need to customize the definition of the (non-copy and non-move) constructors, call
/// `OOPS_CONCRETE_PARAMETERS_ENABLE_COPY_AND_MOVE(className, baseClassName)` instead of this macro
/// and define the constructors on your own.
#define OOPS_CONCRETE_PARAMETERS(className, baseClassName) \
 public: \
  className() : baseClassName() {} \
  explicit className(oops::Parameters* parent) : baseClassName(parent) {} \
  OOPS_CONCRETE_PARAMETERS_ENABLE_COPY_AND_MOVE(className, baseClassName)

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
class Parameters : public ParameterBase, public util::Printable {
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

  void print(std::ostream &os) const override;

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

  /// \brief If OOPS has been built with support for JSON Schema-based validation, check if the
  /// configuration \p config can be deserialized into this Parameters object and if not, throw an
  /// exception with an explanatory message. If validation is not supported, return immediately.
  ///
  /// The check is done by validating the JSON representation of \p config against the JSON schema
  /// returned by jsonSchema().
  ///
  /// To check if JSON Schema-based validation is supported, call isValidationSupported.
  void validate(const eckit::Configuration &config);

  /// \brief Call validate(config) and then deserialize(config).
  void validateAndDeserialize(const eckit::Configuration &config);

  /// \brief Return true if OOPS has been built with support for JSON Schema-based validation (which
  /// requires certain third-party libraries), false otherwise.
  static bool isValidationSupported();

  /// \brief Store the values of all previously registered parameters in \p config.
  void serialize(eckit::LocalConfiguration &config) const override;

  /// \brief Store the values of all previously registered parameters in a new LocalConfiguration
  /// object and return it.
  eckit::LocalConfiguration toConfiguration() const;

  ObjectJsonSchema jsonSchema() const override;

 private:
  std::vector<ParameterBase*> children_;
};

/// \brief Deserialize the configuration \p config into a new instance of \c ParametersType
/// after validating it against the JSON schema defined by \c ParametersType.
///
/// If \p config is not the top-level configuration loaded from a YAML file, then paths to nodes
/// of the YAML tree included in any error messages will be incomplete.
///
/// \tparam ParametersType A concrete subclass of Parameters.
template <typename ParametersType>
ParametersType validateAndDeserialize(const eckit::Configuration &config) {
  ParametersType parameters;
  parameters.validateAndDeserialize(config);
  return parameters;
}

/// \brief Deserialize configurations \p configs into new instances of \c ParametersType
/// after validating them against the JSON schema defined by \c ParametersType.
///
/// If \p configs are not top-level configurations loaded from YAML files, then paths to nodes
/// of the YAML tree included in any error messages will be incomplete.
///
/// \tparam ParametersType A concrete subclass of Parameters.
template <typename ParametersType>
std::vector<ParametersType> validateAndDeserialize(
    const std::vector<eckit::LocalConfiguration> &configs) {
  std::vector<ParametersType> parameters(configs.size());
  for (size_t i = 0; i < configs.size(); ++i)
    parameters[i].validateAndDeserialize(configs[i]);
  return parameters;
}

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERS_H_
