/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_VARIABLECHANGEBASE_H_
#define OOPS_BASE_VARIABLECHANGEBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/make_unique.hpp>
#include <boost/noncopyable.hpp>

#include "oops/base/VariableChangeParametersBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/State.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/PolymorphicParameter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -------------------------------------------------------------------------------------------------

/// Base class for generic variable transform
///
/// Note: subclasses can opt to extract their settings either from a Configuration object or from a
/// subclass of Parameters.
///
/// In the former case, they should provide a constructor taking a const reference to an
/// eckit::Configuration object. In the latter case, the implementer should first define a subclass
/// of Parameters holding the settings of the variable change in question. The latter should
/// then typedef `Parameters_` to the name of that subclass and provide a constructor taking a
/// const reference to an instance of that subclass.
template <typename MODEL>
class VariableChangeBase : public util::Printable,
                           private boost::noncopyable {
  typedef State<MODEL>               State_;

 public:
  explicit VariableChangeBase(const VariableChangeParametersBase &);
  explicit VariableChangeBase(const eckit::Configuration &);
  virtual ~VariableChangeBase() {}

  void setInputVariables(const Variables & vars) { varin_.reset(new Variables(vars)); }
  void setOutputVariables(const Variables & vars) { varout_.reset(new Variables(vars)); }

  virtual void changeVar(const State_ &, State_ &) const = 0;
  virtual void changeVarInverse(const State_ &, State_ &) const = 0;

  State_ changeVar(const State_ &) const;
  State_ changeVarInverse(const State_ &) const;

 private:
  virtual void print(std::ostream &) const = 0;
  std::unique_ptr<Variables> varin_;
  std::unique_ptr<Variables> varout_;
};

// =============================================================================

template <typename MODEL>
class VariableChangeFactory;

// -----------------------------------------------------------------------------

/// \brief A subclass of VariableChangeParametersBase storing the values of all options in a
/// single Configuration object.
///
/// This object can be accessed by calling the value() method of the \p config member variable.
///
/// The ConfigurationParameter class does not perform any parameter validation; models using
/// GenericVariableChangeParameters should therefore ideally be refactored, replacing this
/// class with a dedicated subclass of VariableChangeParametersBase storing each parameter in
/// a separate (Optional/Required)Parameter object.
class GenericVariableChangeParameters : public VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericVariableChangeParameters,
                           VariableChangeParametersBase)
 public:
  ConfigurationParameter config{this};
};

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// VariableChangeParametersBase.
template <typename MODEL>
class VariableChangeParametersWrapper : public Parameters {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of VariableChangeParametersBase
  /// controlling the behavior of a  variable change. The type of the subclass is determined
  /// by the value of the "variable change" key in the Configuration object from which this object
  /// is deserialized.
  PolymorphicParameter<VariableChangeParametersBase, VariableChangeFactory<MODEL>>
    variableChangeParameters{"variable change", "Identity", this};
};

// =============================================================================

/// VariableChange factory
template <typename MODEL>
class VariableChangeFactory {
  typedef Geometry<MODEL>   Geometry_;

 public:
  /// \brief Create and return a new variable change.
  ///
  /// The type of the variable change is determined by the `variable change` attribute of \p
  /// parameters. \p parameters must be an instance of the subclass of
  /// VariableChangeParametersBase associated with that variable change type, otherwise an
  /// exception will be thrown.
  static VariableChangeBase<MODEL> * create(const VariableChangeParametersBase &,
                                            const Geometry_ &);

  /// \brief Create and return a new variable change.
  ///
  /// Deprecated overload taking a Configuration instead of a VariableChangeParametersBase.
  static VariableChangeBase<MODEL> * create(const eckit::Configuration &, const Geometry_ &);

  /// \brief Create and return an instance of the subclass of VariableChangeParametersBase
  /// storing parameters of variable changes of the specified type.
  static std::unique_ptr<VariableChangeParametersBase> createParameters(const std::string &name);

  /// \brief Return the names of all variable changes that can be created by one of the
  /// registered makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~VariableChangeFactory() = default;

 protected:
  /// \brief Register a maker able to create variable changes of type \p name.
  explicit VariableChangeFactory(const std::string &);

 private:
  virtual VariableChangeBase<MODEL> * make(const VariableChangeParametersBase &,
                                           const Geometry_ &) = 0;

  virtual std::unique_ptr<VariableChangeParametersBase> makeParameters() const = 0;

  static std::map < std::string, VariableChangeFactory<MODEL> * > & getMakers() {
    static std::map < std::string, VariableChangeFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -------------------------------------------------------------------------------------------------

template<class MODEL, class T>
class VariableChangeMaker : public VariableChangeFactory<MODEL> {
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// GenericVariableChangeParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<T, GenericVariableChangeParameters> Parameters_;

  typedef Geometry<MODEL>   Geometry_;

  VariableChangeBase<MODEL> * make(const VariableChangeParametersBase & params,
                                   const Geometry_ & resol) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(resol,
                 parametersOrConfiguration<HasParameters_<T>::value>(stronglyTypedParams));
  }

  std::unique_ptr<VariableChangeParametersBase> makeParameters() const override {
    return boost::make_unique<Parameters_>();
  }

 public:
  explicit VariableChangeMaker(const std::string & name)
    : VariableChangeFactory<MODEL>(name) {}
};

// -------------------------------------------------------------------------------------------------

template <typename MODEL>
VariableChangeFactory<MODEL>::VariableChangeFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in the variable change factory.");
  }
  getMakers()[name] = this;
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL>
VariableChangeBase<MODEL> * VariableChangeFactory<MODEL>::create(
    const VariableChangeParametersBase & params, const Geometry_ & resol)
{
  Log::trace() << "VariableChangeBase<MODEL>::create starting" << std::endl;
// Not good: should not create anything if no variable change required. YT
  const std::string &id = params.variableChange.value().value();
  typename std::map<std::string, VariableChangeFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in variable change factory.");
  }
  VariableChangeBase<MODEL> * ptr = jerr->second->make(params, resol);
  Log::trace() << "VariableChangeBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -------------------------------------------------------------------------------------------------

template <typename MODEL>
VariableChangeBase<MODEL> * VariableChangeFactory<MODEL>::create(const eckit::Configuration & conf,
                                                                 const Geometry_ & resol) {
  VariableChangeParametersWrapper<MODEL> parameters;
  parameters.validateAndDeserialize(conf);
  return create(parameters.variableChangeParameters, resol);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<VariableChangeParametersBase> VariableChangeFactory<MODEL>::createParameters(
    const std::string &name) {
  typename std::map<std::string, VariableChangeFactory<MODEL>*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in VariableChangeFactory");
  }
  return it->second->makeParameters();
}

// =================================================================================================

template<typename MODEL>
VariableChangeBase<MODEL>::VariableChangeBase(const VariableChangeParametersBase & params)
  : varin_(), varout_()
{
  if (params.inputVariables.value() != boost::none) {
    varin_.reset(new Variables(*params.inputVariables.value()));
    Log::trace() << "VariableChangeBase<MODEL>::VariableChangeBase input variables: "
                 << *varin_ << std::endl;
  }
  if (params.outputVariables.value() != boost::none) {
    varout_.reset(new Variables(*params.outputVariables.value()));
    Log::trace() << "VariableChangeBase<MODEL>::VariableChangeBase output variables: "
                 << *varout_ << std::endl;
  }
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
VariableChangeBase<MODEL>::VariableChangeBase(const eckit::Configuration & conf)
  : VariableChangeBase(validateAndDeserialize<GenericVariableChangeParameters>(conf))
{}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> VariableChangeBase<MODEL>::changeVar(const State_ & xin) const {
  Log::trace() << "VariableChangeBase<MODEL>::changeVar starting" << std::endl;
  ASSERT(varout_);
  State_ xout(xin.geometry(), *varout_, xin.validTime());
  this->changeVar(xin, xout);
  Log::trace() << "VariableChangeBase<MODEL>::changeVar done" << std::endl;
  return xout;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> VariableChangeBase<MODEL>::changeVarInverse(const State_ & xin) const {
  Log::trace() << "VariableChangeBase<MODEL>::changeVarInverse starting" << std::endl;
  ASSERT(varin_);
  State_ xout(xin.geometry(), *varin_, xin.validTime());
  this->changeVarInverse(xin, xout);
  Log::trace() << "VariableChangeBase<MODEL>::changeVarInverse done" << std::endl;
  return xout;
}

// -------------------------------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_VARIABLECHANGEBASE_H_
