/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_LINEARVARIABLECHANGEBASE_H_
#define OOPS_BASE_LINEARVARIABLECHANGEBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/make_unique.hpp>
#include <boost/noncopyable.hpp>

#include "eckit/exception/Exceptions.h"

#include "oops/base/LinearVariableChangeParametersBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// Base class for generic variable transform
///
/// Note: subclasses can opt to extract their settings either from a Configuration object or from a
/// subclass of Parameters.
///
/// In the former case, they should provide a constructor taking a const reference to an
/// eckit::Configuration object. In the latter case, the implementer should first define a subclass
/// of Parameters holding the settings of the linear variable change in question. The latter should
/// then typedef `Parameters_` to the name of that subclass and provide a constructor taking a
/// const reference to an instance of that subclass.
template <typename MODEL>
class LinearVariableChangeBase : public util::Printable,
                                 private boost::noncopyable {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  explicit LinearVariableChangeBase(const LinearVariableChangeParametersBase &);
  explicit LinearVariableChangeBase(const eckit::Configuration &);
  virtual ~LinearVariableChangeBase() {}

  void setInputVariables(const Variables & vars) { varin_.reset(new Variables(vars)); }
  void setOutputVariables(const Variables & vars) { varout_.reset(new Variables(vars)); }

  virtual void multiply(const Increment_ &, Increment_ &) const = 0;
  virtual void multiplyInverse(const Increment_ &, Increment_ &) const = 0;
  virtual void multiplyAD(const Increment_ &, Increment_ &) const = 0;
  virtual void multiplyInverseAD(const Increment_ &, Increment_ &) const = 0;

  Increment_ multiply(const Increment_ &) const;
  Increment_ multiplyInverse(const Increment_ &) const;
  Increment_ multiplyAD(const Increment_ &) const;
  Increment_ multiplyInverseAD(const Increment_ &) const;

 private:
  virtual void print(std::ostream &) const = 0;
  std::unique_ptr<Variables> varin_;
  std::unique_ptr<Variables> varout_;
};

// =============================================================================

template <typename MODEL>
class LinearVariableChangeFactory;

// -----------------------------------------------------------------------------

/// \brief A subclass of LinearVariableChangeParametersBase storing the values of all options in a
/// single Configuration object.
///
/// This object can be accessed by calling the value() method of the \p config member variable.
///
/// The ConfigurationParameter class does not perform any parameter validation; models using
/// GenericLinearVariableChangeParameters should therefore ideally be refactored, replacing this
/// class with a dedicated subclass of LinearVariableChangeParametersBase storing each parameter in
/// a separate (Optional/Required)Parameter object.
class GenericLinearVariableChangeParameters : public LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericLinearVariableChangeParameters,
                           LinearVariableChangeParametersBase)
 public:
  ConfigurationParameter config{this};
};

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// LinearVariableChangeParametersBase.
template <typename MODEL>
class LinearVariableChangeParametersWrapper : public Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of LinearVariableChangeParametersBase
  /// controlling the behavior of a linear variable change. The type of the subclass is determined
  /// by the value of the "variable change" key in the Configuration object from which this object
  /// is deserialized.
  RequiredPolymorphicParameter<LinearVariableChangeParametersBase,
                               LinearVariableChangeFactory<MODEL>>
    variableChangeParameters{"variable change", this};
};

// =============================================================================

/// LinearVariableChange factory
template <typename MODEL>
class LinearVariableChangeFactory {
  typedef Geometry<MODEL>   Geometry_;
  typedef State<MODEL>      State_;

 public:
  /// \brief Create and return a new linear variable change.
  ///
  /// The type of the variable change is determined by the `variable change` attribute of \p
  /// parameters. \p parameters must be an instance of the subclass of
  /// LinearVariableChangeParametersBase associated with that variable change type, otherwise an
  /// exception will be thrown.
  static LinearVariableChangeBase<MODEL> * create(
      const State_ &, const State_ &, const Geometry_ &,
      const LinearVariableChangeParametersBase & parameters);

  /// \brief Create and return a new linear variable change.
  ///
  /// Deprecated overload taking a Configuration instead of a LinearVariableChangeParametersBase.
  static LinearVariableChangeBase<MODEL> * create(
      const State_ &, const State_ &, const Geometry_ &,
      const eckit::Configuration & conf);

  /// \brief Create and return an instance of the subclass of LinearVariableChangeParametersBase
  /// storing parameters of linear variable changes of the specified type.
  static std::unique_ptr<LinearVariableChangeParametersBase> createParameters(
      const std::string &name);

  /// \brief Return the names of all linear variable changes that can be created by one of the
  /// registered makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~LinearVariableChangeFactory() = default;

 protected:
  /// \brief Register a maker able to create linear variable changes of type \p name.
  explicit LinearVariableChangeFactory(const std::string &name);

 private:
  virtual LinearVariableChangeBase<MODEL> * make(const State_ &, const State_ &,
                                                 const Geometry_ &,
                                                 const LinearVariableChangeParametersBase &) = 0;

  virtual std::unique_ptr<LinearVariableChangeParametersBase> makeParameters() const = 0;

  static std::map < std::string, LinearVariableChangeFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LinearVariableChangeFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LinearVariableChangeMaker : public LinearVariableChangeFactory<MODEL> {
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// GenericLinearVariableChangeParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<T, GenericLinearVariableChangeParameters>
    Parameters_;

  typedef Geometry<MODEL>   Geometry_;
  typedef State<MODEL>      State_;

  LinearVariableChangeBase<MODEL> * make(const State_ & bg, const State_ & fg,
                                       const Geometry_ & geom,
                                       const LinearVariableChangeParametersBase& params) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(bg, fg, geom,
                 parametersOrConfiguration<HasParameters_<T>::value>(stronglyTypedParams));
  }

  std::unique_ptr<LinearVariableChangeParametersBase> makeParameters() const override {
    return boost::make_unique<Parameters_>();
  }

 public:
  explicit LinearVariableChangeMaker(const std::string & name)
    : LinearVariableChangeFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearVariableChangeFactory<MODEL>::LinearVariableChangeFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in the linear variable change factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearVariableChangeBase<MODEL> * LinearVariableChangeFactory<MODEL>::create(
     const State_ & bg, const State_ & fg,
     const Geometry_ & geom, const LinearVariableChangeParametersBase & parameters) {
  Log::trace() << "LinearVariableChangeBase<MODEL>::create starting" << std::endl;
  const std::string &id = parameters.variableChange.value().value();
  typename std::map<std::string, LinearVariableChangeFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);

  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in LinearVariableChangeFactory." << std::endl;
    Log::error() << "Factory contains " << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, LinearVariableChangeFactory<MODEL>*>::const_iterator
         jj = getMakers().begin(); jj !=  getMakers().end(); ++jj) {
      Log::error() << "A " << jj->first << " variable change option" << std::endl;
    }
    throw std::runtime_error(id + " does not exist in the linear variable change factory.");
  }

  LinearVariableChangeBase<MODEL> * ptr = jerr->second->make(bg, fg, geom, parameters);
  Log::trace() << "LinearVariableChangeBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearVariableChangeBase<MODEL> * LinearVariableChangeFactory<MODEL>::create(
     const State_ & bg, const State_ & fg,
     const Geometry_ & geom, const eckit::Configuration & conf) {
  LinearVariableChangeParametersWrapper<MODEL> parameters;
  parameters.validateAndDeserialize(conf);
  return create(bg, fg, geom, parameters.variableChangeParameters);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<LinearVariableChangeParametersBase>
  LinearVariableChangeFactory<MODEL>::createParameters(
    const std::string &name) {
  typename std::map<std::string, LinearVariableChangeFactory<MODEL>*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in LinearVariableChangeFactory");
  }
  return it->second->makeParameters();
}

// =============================================================================

template<typename MODEL>
LinearVariableChangeBase<MODEL>::LinearVariableChangeBase(
    const LinearVariableChangeParametersBase & parameters)
  : varin_(), varout_()
{
  if (parameters.inputVariables.value() != boost::none) {
    varin_.reset(new Variables(*parameters.inputVariables.value()));
    Log::trace() << "LinearVariableChangeBase<MODEL>::LinearVariableChangeBase input variables: "
                 << *varin_ << std::endl;
  }
  if (parameters.outputVariables.value() != boost::none) {
    varout_.reset(new Variables(*parameters.outputVariables.value()));
    Log::trace() << "LinearVariableChangeBase<MODEL>::LinearVariableChangeBase output variables: "
                 << *varout_ << std::endl;
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LinearVariableChangeBase<MODEL>::LinearVariableChangeBase(
    const eckit::Configuration & conf)
  : LinearVariableChangeBase(validateAndDeserialize<GenericLinearVariableChangeParameters>(conf))
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> LinearVariableChangeBase<MODEL>::multiply(const Increment<MODEL> & dxin) const {
  ASSERT(varin_);
  Increment_ dxout(dxin.geometry(), *varout_, dxin.validTime());
  this->multiply(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> LinearVariableChangeBase<MODEL>::multiplyAD(const Increment_ & dxin) const {
  ASSERT(varout_);
  Increment_ dxout(dxin.geometry(), *varin_, dxin.validTime());
  this->multiplyAD(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> LinearVariableChangeBase<MODEL>::multiplyInverse(const Increment_ & dxin) const {
  ASSERT(varout_);
  Increment_ dxout(dxin.geometry(), *varin_, dxin.validTime());
  this->multiplyInverse(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> LinearVariableChangeBase<MODEL>::multiplyInverseAD(const Increment_ & dxin) const {
  ASSERT(varin_);
  Increment_ dxout(dxin.geometry(), *varout_, dxin.validTime());
  this->multiplyInverseAD(dxin, dxout);
  return dxout;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_LINEARVARIABLECHANGEBASE_H_
