/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERSORCONFIGURATION_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERSORCONFIGURATION_H_

namespace eckit {
class Configuration;
class LocalConfiguration;
}

namespace oops {

/// \brief Used to implement the parametersOrConfiguration() function. Shouldn't need to be used
/// elsewhere.
///
/// It's needed because C++ doesn't support partial specialization for function templates, only
/// for class templates.
template <bool returnParameters, typename ParametersType>
struct ParametersOrConfiguration;

template <typename ParametersType>
struct ParametersOrConfiguration<true /* returnParameters? */, ParametersType> {
  typedef ParametersType ReturnType;
  const ReturnType & operator()(const ParametersType& parameters) const {
    return parameters;
  }
};

template <typename ParametersType>
struct ParametersOrConfiguration<false /* returnParameters? */, ParametersType> {
  typedef eckit::LocalConfiguration ReturnType;
  const ReturnType & operator()(const ParametersType& parameters) const {
    return parameters.config.value();
  }
};

/// \brief Helper function used to pass arguments to functions that may take a const reference to
/// either eckit::Configuration ("old style") or to a subclass to Parameters ("new style").
///
/// \tparam returnParameters
///   If true, this function will simply return its argument. If false, it will assume the argument
///   is a subclass of Parameters containing a `config` member of type ConfigurationParameter, and
///   it will return the reference to the LocalConfiguration object stored in that parameter.
/// \tparam ParametersType
///   A subclass of Parameters.
template <bool returnParameters, typename ParametersType>
const typename ParametersOrConfiguration<returnParameters, ParametersType>::ReturnType &
parametersOrConfiguration(const ParametersType& parameters) {
  return ParametersOrConfiguration<returnParameters, ParametersType>()(parameters);
}

/// Overloaded function used to handle in a uniform way "new-style" interface implementations that
/// provide functions returning instances of subclasses of Parameters and "old style"
/// implementations that provide functions returning eckit::Configuration objects.
///
/// This is the overload handling the former of these cases; it simply returns its argument.
template <typename ParametersType>
const ParametersType& toParameters(const ParametersType &parameters) {
  return parameters;
}

/// Overloaded function used to handle in a uniform way interface implementations that
/// provide functions returning instances of subclasses of Parameters and those that provide
/// functions returning eckit::Configuration objects.
///
/// This is the overload handling the latter of these cases; it deserializes the received
/// Configuration object into an instance of ParametersType (which would typically be
/// GenericParameters or a similar class) and returns that instance.
template <typename ParametersType>
ParametersType toParameters(const eckit::Configuration &conf) {
  ParametersType parameters;
  parameters.validateAndDeserialize(conf);
  return parameters;
}

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERSORCONFIGURATION_H_
