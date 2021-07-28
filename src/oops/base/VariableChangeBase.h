/*
 * (C) Copyright 2018-2021 UCAR
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

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/VariableChangeParametersBase.h"
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

/// Base class for a variable transforms, defining the interfaces.
/// Use this class as a base class for generic implementations,
/// and VariableChangeBase as a base class for MODEL-specific implementations.
///
/// Note: subclasses can opt to extract their settings either from a Configuration object or from a
/// subclass of Parameters.
///
/// In the former case, they should provide a constructor with the following signature:
///
///    VariableChange(const Geometry_ &, const eckit::Configuration &);
///
/// In the latter case, the implementer should first define a subclass of
/// VariableChangeParametersBase holding the settings of the variable change in question.
/// The implementation of the VariableChange interface should then typedef `Parameters_`
/// to the name of that subclass and provide a constructor with the following signature:
///
///    VariableChange(const Geometry_ &, const Parameters_ &);
template <typename MODEL>
class GenericVariableChangeBase : public util::Printable,
                                  private boost::noncopyable {
  typedef State<MODEL>    State_;

 public:
  GenericVariableChangeBase() = default;
  virtual ~GenericVariableChangeBase() = default;

  /// change variables from state \p xin to \p xout
  virtual void changeVar(const State_ & xin, State_ & xout) const = 0;
  /// inverse of changeVar, change variables back from \p xout to \p xin
  virtual void changeVarInverse(const State_ & xout, State_ & xin) const = 0;

 private:
  /// Print, used for logging
  virtual void print(std::ostream &) const = 0;
};

/// \brief Base class for MODEL-specific implementations of VariableChange class.
/// The complete interface that needs to be implemented is described in
/// GenericVariableChangeBase. VariableChangeBase overrides GenericVariableChangeBase
/// methods to pass MODEL-specific implementations of State to the MODEL-specific
/// implementation of VariableChange.
template <typename MODEL>
class VariableChangeBase : public GenericVariableChangeBase<MODEL> {
  typedef typename MODEL::State    State_;

 public:
  VariableChangeBase() = default;
  virtual ~VariableChangeBase() = default;

  /// Overrides for VariableChangeBase classes, passing MODEL-specific classes to the
  /// MODEL-specific implementations of VariableChange
  void changeVar(const State<MODEL> & xin, State<MODEL> & xout) const final
    { this->changeVar(xin.state(), xout.state()); }
  void changeVarInverse(const State<MODEL> & xout, State<MODEL> & xin) const final
    { this->changeVarInverse(xout.state(), xin.state()); }

  /// change variables from state \p xin to \p xout
  virtual void changeVar(const State_ & xin, State_ & xout) const = 0;
  /// inverse of changeVar, change variables back from \p xout to \p xin
  virtual void changeVarInverse(const State_ & xout, State_ & xin) const = 0;
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
  static GenericVariableChangeBase<MODEL> * create(const Geometry_ &,
                                                   const VariableChangeParametersBase &);
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
  virtual GenericVariableChangeBase<MODEL> * make(const Geometry_ &,
                                                  const VariableChangeParametersBase &) = 0;

  virtual std::unique_ptr<VariableChangeParametersBase> makeParameters() const = 0;

  static std::map < std::string, VariableChangeFactory<MODEL> * > & getMakers() {
    static std::map < std::string, VariableChangeFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -------------------------------------------------------------------------------------------------

template<class MODEL, class T>
class GenericVariableChangeMaker : public VariableChangeFactory<MODEL> {
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// GenericVariableChangeParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<T, GenericVariableChangeParameters> Parameters_;

  typedef Geometry<MODEL>   Geometry_;

  GenericVariableChangeBase<MODEL> * make(const Geometry_ & resol,
                                          const VariableChangeParametersBase & params) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(resol,
                 parametersOrConfiguration<HasParameters_<T>::value>(stronglyTypedParams));
  }

  std::unique_ptr<VariableChangeParametersBase> makeParameters() const override {
    return boost::make_unique<Parameters_>();
  }

 public:
  explicit GenericVariableChangeMaker(const std::string & name)
    : VariableChangeFactory<MODEL>(name) {}
};


// -------------------------------------------------------------------------------------------------

template<class MODEL, class T>
class VariableChangeMaker : public VariableChangeFactory<MODEL> {
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// GenericVariableChangeParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<T, GenericVariableChangeParameters> Parameters_;

  typedef Geometry<MODEL>   Geometry_;

  VariableChangeBase<MODEL> * make(const Geometry_ & resol,
                                   const VariableChangeParametersBase & params) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(resol.geometry(),
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
GenericVariableChangeBase<MODEL> * VariableChangeFactory<MODEL>::create(
    const Geometry_ & resol, const VariableChangeParametersBase & params)
{
  Log::trace() << "VariableChangeBase<MODEL>::create starting" << std::endl;
// Not good: should not create anything if no variable change required. YT
  const std::string &id = params.variableChange.value();
  typename std::map<std::string, VariableChangeFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in variable change factory.");
  }
  GenericVariableChangeBase<MODEL> * ptr = jerr->second->make(resol, params);
  Log::trace() << "VariableChangeBase<MODEL>::create done" << std::endl;
  return ptr;
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

// -------------------------------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_VARIABLECHANGEBASE_H_
