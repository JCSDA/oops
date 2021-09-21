/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_LINEARMODELBASE_H_
#define OOPS_GENERIC_LINEARMODELBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/make_unique.hpp>
#include <boost/noncopyable.hpp>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Base class for generic implementations of the linearized forecasting models.
/// Use this class as a base class for generic implementations,
/// and interface::LinearModelBase as a base class for MODEL-specific implementations.
///
/// Note: implementations of this interface can opt to extract their settings either from
/// a Configuration object or from a subclass of LinearModelParametersBase.
///
/// In the former case, they should provide a constructor with the following signature:
///
///    LinearModelBase(const Geometry_ &, const eckit::Configuration &);
///
/// In the latter case, the implementer should first define a subclass of LinearModelParametersBase
/// holding the settings of the linear model in question. The implementation of the LinearModelBase
/// interface should then typedef `Parameters_` to the name of that subclass and provide a
/// constructor with the following signature:
///
///    LinearModelBase(const Geometry_ &, const Parameters_ &);
///
template <typename MODEL>
class LinearModelBase : public util::Printable,
                        private boost::noncopyable {
  typedef Increment<MODEL>         Increment_;
  typedef ModelAuxControl<MODEL>   ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL> ModelAuxInc_;
  typedef State<MODEL>             State_;

 public:
  static const std::string classname() {return "oops::LinearModelBase";}

  LinearModelBase() = default;
  virtual ~LinearModelBase() = default;

  /// \brief Tangent linear initialization, called before every run
  virtual void initializeTL(Increment_ &) const = 0;
  /// \brief Tangent linear "step", called during run; updates increment to the next time
  virtual void stepTL(Increment_ &, const ModelAuxInc_ &) const = 0;
  /// \brief Tangent linear finalization; called after each run
  virtual void finalizeTL(Increment_ &) const = 0;

  /// \brief Tangent linear initialization, called before every run
  virtual void initializeAD(Increment_ &) const = 0;
  /// \brief Tangent linear "step", called during run; updates increment to the next time
  virtual void stepAD(Increment_ &, ModelAuxInc_ &) const = 0;
  /// \brief Tangent linear finalization; called after each run
  virtual void finalizeAD(Increment_ &) const = 0;

  /// \brief Set the trajectory for the linear model, called after each step of the forecast
  virtual void setTrajectory(const State_ &, State_ &, const ModelAuxCtl_ &) = 0;

  /// \brief Time step for running LinearModel's forecast in oops (frequency with which the
  /// increment will be updated)
  virtual const util::Duration & timeResolution() const = 0;
  /// \brief LinearModel variables (only used in 4DVar)
  virtual const oops::Variables & variables() const = 0;

 private:
  /// \brief Print; used for logging
  virtual void print(std::ostream &) const = 0;
};

// =============================================================================

template <typename MODEL>
class LinearModelFactory;

// -----------------------------------------------------------------------------

/// \brief Base class for classes storing linear model-specific parameters.
class LinearModelParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(LinearModelParametersBase, Parameters)
 public:
  /// \brief LinearModel name.
  ///
  /// \note This parameter is marked as optional because it is only required in certain
  /// circumstances (e.g. when linear model parameters are deserialized into a
  /// LinearModelParametersWrapper and used by LinearModelFactory to instantiate a linear model
  /// whose type is determined at runtime), but not others (e.g. in tests written with a particular
  /// linear model in mind). LinearModelParametersWrapper will throw an exception if this parameter
  /// is not provided.
  OptionalParameter<std::string> name{"name", this};
};

// -----------------------------------------------------------------------------

/// \brief A subclass of LinearModelParametersBase storing the values of all options in a
/// single Configuration object.
///
/// This object can be accessed by calling the value() method of the \p config member variable.
///
/// The ConfigurationParameter class does not perform any parameter validation; linear models using
/// GenericLinearModelParameters should therefore ideally be refactored, replacing this class with a
/// dedicated subclass of LinearModelParametersBase storing each parameter in a separate
/// (Optional/Required)Parameter object.
class GenericLinearModelParameters : public LinearModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericLinearModelParameters, LinearModelParametersBase)
 public:
  ConfigurationParameter config{this};
};

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// LinearModelParametersBase.
template <typename MODEL>
class LinearModelParametersWrapper : public Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearModelParametersWrapper, Parameters)
 public:
  /// After deserializtion, holds an instance of a subclass of LinearModelParametersBase controlling
  /// the behavior of a linear model. The type of the subclass is determined by the value of the
  /// "name" key in the Configuration object from which this object is deserialized.
  RequiredPolymorphicParameter<LinearModelParametersBase, LinearModelFactory<MODEL>>
    linearModelParameters{"name", this};
};

// =============================================================================

/// LinearModel factory
template <typename MODEL>
class LinearModelFactory {
  typedef Geometry<MODEL>   Geometry_;

 public:
  /// \brief Create and return a new linear model.
  ///
  /// The linear model's type is determined by the \c name attribute of \p parameters.
  /// \p parameters must be an instance of the subclass of LinearModelParametersBase
  /// associated with that linear model type, otherwise an exception will be thrown.
  static LinearModelBase<MODEL> * create(const Geometry_ &,
                                         const LinearModelParametersBase & parameters);

  /// \brief Create and return an instance of the subclass of LinearModelParametersBase
  /// storing parameters of linear models of the specified type.
  static std::unique_ptr<LinearModelParametersBase> createParameters(const std::string &name);

  /// \brief Return the names of all linear models that can be created by one of the registered
  /// makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~LinearModelFactory() = default;

 protected:
  /// \brief Register a maker able to create linear models of type \p name.
  explicit LinearModelFactory(const std::string & name);

 private:
  virtual LinearModelBase<MODEL> * make(const Geometry_ &, const LinearModelParametersBase &) = 0;

  virtual std::unique_ptr<LinearModelParametersBase> makeParameters() const = 0;

  static std::map < std::string, LinearModelFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LinearModelFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

/// \brief A subclass of LinearModelFactory able to create instances of T (a concrete subclass of
/// LinearModelBase<MODEL>). Passes Geometry<MODEL> to the constructor of T.
template<class MODEL, class T>
class LinearModelMaker : public LinearModelFactory<MODEL> {
 private:
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// GenericLinearModelParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<T, GenericLinearModelParameters> Parameters_;

 public:
  typedef Geometry<MODEL>   Geometry_;

  explicit LinearModelMaker(const std::string & name) : LinearModelFactory<MODEL>(name) {}

  LinearModelBase<MODEL> * make(const Geometry_ & geom,
                                const LinearModelParametersBase & parameters) override {
    Log::trace() << "LinearModelBase<MODEL>::make starting" << std::endl;
    const auto &stronglyTypedParameters = dynamic_cast<const Parameters_&>(parameters);
    return new T(geom,
                 parametersOrConfiguration<HasParameters_<T>::value>(stronglyTypedParameters));
  }

  std::unique_ptr<LinearModelParametersBase> makeParameters() const override {
    return boost::make_unique<Parameters_>();
  }
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearModelFactory<MODEL>::LinearModelFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in the linear model factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearModelBase<MODEL> * LinearModelFactory<MODEL>::create(const Geometry_ & geom,
                                                    const LinearModelParametersBase & parameters) {
  Log::trace() << "LinearModelFactory<MODEL>::create starting" << std::endl;
  const std::string &id = parameters.name.value().value();
  typename std::map<std::string, LinearModelFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in the linear model factory");
  }
  LinearModelBase<MODEL> * ptr = jerr->second->make(geom, parameters);
  Log::trace() << "LinearModelFactory<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<LinearModelParametersBase> LinearModelFactory<MODEL>::createParameters(
    const std::string &name) {
  Log::trace() << "LinearModelFactory<MODEL>::createParameters starting" << std::endl;
  typename std::map<std::string, LinearModelFactory<MODEL>*>::iterator it = getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in the linear model factory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_LINEARMODELBASE_H_
