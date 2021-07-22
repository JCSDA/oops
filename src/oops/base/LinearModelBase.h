/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_LINEARMODELBASE_H_
#define OOPS_BASE_LINEARMODELBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/make_unique.hpp>
#include <boost/noncopyable.hpp>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/interface/State.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// Base class for encapsulation of the linear forecast model.
/*!
 * Defines the interfaces for the linear model.
 */

// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearModelBase : public util::Printable,
                        private boost::noncopyable {
  typedef Increment<MODEL>              Increment_;
  typedef Geometry<MODEL>               Geometry_;
  typedef ModelAuxControl<MODEL>        ModelAux_;
  typedef ModelAuxIncrement<MODEL>      ModelAuxIncr_;
  typedef State<MODEL>                  State_;

 public:
  static const std::string classname() {return "oops::LinearModelBase";}

  LinearModelBase() {}
  virtual ~LinearModelBase() {}

// Set the linearization trajectory
  void setTrajectory(const State_ &, State_ &, const ModelAux_ &);

// Run the TL forecast
  void initializeTL(Increment_ &) const;
  void stepTL(Increment_ &, const ModelAuxIncr_ &) const;
  void finalizeTL(Increment_ &) const;

// Run the AD forecast
  void initializeAD(Increment_ &) const;
  void stepAD(Increment_ &, ModelAuxIncr_ &) const;
  void finalizeAD(Increment_ &) const;

// Information and diagnostics
  virtual const util::Duration & timeResolution() const = 0;
  virtual const oops::Variables & variables() const = 0;

 protected:
// Set the linearization trajectory
  virtual void setTrajectory(const typename MODEL::State &, typename MODEL::State &,
                             const typename MODEL::ModelAuxControl &) = 0;

// Run the TL forecast
  virtual void initializeTL(typename MODEL::Increment &) const = 0;
  virtual void stepTL(typename MODEL::Increment &,
                      const typename MODEL::ModelAuxIncrement &) const = 0;
  virtual void finalizeTL(typename MODEL::Increment &) const = 0;

// Run the AD forecast
  virtual void initializeAD(typename MODEL::Increment &) const = 0;
  virtual void stepAD(typename MODEL::Increment &, typename MODEL::ModelAuxIncrement &) const = 0;
  virtual void finalizeAD(typename MODEL::Increment &) const = 0;

// Information and diagnostics
  virtual void print(std::ostream &) const = 0;
};

// =============================================================================

template <typename MODEL>
class LinearModelFactory;

// -----------------------------------------------------------------------------

/// \brief Base class for classes storing model-specific parameters.
class LinearModelParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(LinearModelParametersBase, Parameters)
 public:
  /// \brief Model name.
  ///
  /// \note This parameter is marked as optional because it is only required in certain
  /// circumstances (e.g. when model parameters are deserialized into a
  /// LinearModelParametersWrapper and used by LinearModelFactory to instantiate a tangent linear
  /// model whose type is determined at runtime), but not others (e.g. in tests written with a
  /// particular model in mind). LinearModelParametersWrapper will throw an exception if this
  /// parameter is not provided.
  OptionalParameter<std::string> name{"name", this};
};

// -----------------------------------------------------------------------------

/// \brief A subclass of LinearModelParametersBase storing the values of all options in a
/// single Configuration object.
///
/// This object can be accessed by calling the value() method of the \p config member variable.
///
/// The ConfigurationParameter class does not perform any parameter validation; models using
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
  /// After deserialization, holds an instance of a subclass of LinearModelParametersBase
  /// controlling the behavior of a tangent linear model. The type of the subclass is determined by
  /// the value of the "name" key in the Configuration object from which this object is
  /// deserialized.
  RequiredPolymorphicParameter<LinearModelParametersBase, LinearModelFactory<MODEL>>
    modelParameters{"name", this};
};

// =============================================================================

/// \brief Tangent linear model factory.
template <typename MODEL>
class LinearModelFactory {
  typedef Geometry<MODEL>   Geometry_;

 public:
  /// \brief Create and return a new tangent linear model.
  ///
  /// The model's type is determined by the \c name attribute of \p parameters.
  /// \p parameters must be an instance of the subclass of LinearModelParametersBase
  /// associated with that model type, otherwise an exception will be thrown.
  static LinearModelBase<MODEL> * create(const Geometry_ &,
                                         const LinearModelParametersBase & parameters);

  /// \brief Create and return an instance of the subclass of LinearModelParametersBase
  /// storing parameters of tangent linear models of the specified type.
  static std::unique_ptr<LinearModelParametersBase> createParameters(const std::string &name);

  /// \brief Return the names of all tangent linear models that can be created by one of the
  /// registered makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~LinearModelFactory() = default;

 protected:
  /// \brief Register a maker able to create tangent linear models of type \p name.
  explicit LinearModelFactory(const std::string &name);

 private:
  virtual LinearModelBase<MODEL> * make(const Geometry_ &, const LinearModelParametersBase &) = 0;

  virtual std::unique_ptr<LinearModelParametersBase> makeParameters() const = 0;

  static std::map < std::string, LinearModelFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LinearModelFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LinearModelMaker : public LinearModelFactory<MODEL> {
 private:
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// GenericLinearModelParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<T, GenericLinearModelParameters> Parameters_;

  typedef Geometry<MODEL>   Geometry_;

  LinearModelBase<MODEL> * make(const Geometry_ & geom,
                                const LinearModelParametersBase & parameters) override {
    const auto &stronglyTypedParameters = dynamic_cast<const Parameters_&>(parameters);
    return new T(geom.geometry(),
                 parametersOrConfiguration<HasParameters_<T>::value>(stronglyTypedParameters));
  }

  std::unique_ptr<LinearModelParametersBase> makeParameters() const override {
    return boost::make_unique<Parameters_>();
  }

 public:
  explicit LinearModelMaker(const std::string & name) : LinearModelFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearModelFactory<MODEL>::LinearModelFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in tangent linear model factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearModelBase<MODEL>* LinearModelFactory<MODEL>::create(
    const Geometry_ & geom, const LinearModelParametersBase & parameters) {
  Log::trace() << "LinearModelBase<MODEL>::create starting" << std::endl;
  const std::string &id = parameters.name.value().value();
  typename std::map<std::string, LinearModelFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in the tangent linear model factory." << std::endl;
    Log::error() << "Factory contains " << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, LinearModelFactory<MODEL>*>::const_iterator
         jj = getMakers().begin(); jj !=  getMakers().end(); ++jj) {
      Log::error() << "A " << jj->first << " linear model" << std::endl;
    }
    throw std::runtime_error(id + " does not exist in tangent linear model factory.");
  }
  LinearModelBase<MODEL> * ptr = jerr->second->make(geom, parameters);
  Log::trace() << "LinearModelBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<LinearModelParametersBase> LinearModelFactory<MODEL>::createParameters(
    const std::string &name) {
  typename std::map<std::string, LinearModelFactory<MODEL>*>::iterator it = getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in the tangent linear model factory");
  }
  return it->second->makeParameters();
}

// =============================================================================

template<typename MODEL>
void LinearModelBase<MODEL>::setTrajectory(const State_ & xx, State_ & xlr,
                                          const ModelAux_ & maux) {
  Log::trace() << "LinearModelBase<MODEL>::setTrajectory starting" << std::endl;
  this->setTrajectory(xx.state(), xlr.state(), maux.modelauxcontrol());
  Log::trace() << "LinearModelBase<MODEL>::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::initializeTL(Increment_ & dx) const {
  Log::trace() << "LinearModelBase<MODEL>::initializeTL starting" << std::endl;
  this->initializeTL(dx.increment());
  Log::trace() << "LinearModelBase<MODEL>::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::stepTL(Increment_ & dx, const ModelAuxIncr_ & merr) const {
  Log::trace() << "LinearModelBase<MODEL>::stepTL starting" << std::endl;
  this->stepTL(dx.increment(), merr.modelauxincrement());
  Log::trace() << "LinearModelBase<MODEL>::stepTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::finalizeTL(Increment_ & dx) const {
  Log::trace() << "LinearModelBase<MODEL>::finalizeTL starting" << std::endl;
  this->finalizeTL(dx.increment());
  Log::trace() << "LinearModelBase<MODEL>::finalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::initializeAD(Increment_ & dx) const {
  Log::trace() << "LinearModelBase<MODEL>::initializeAD starting" << std::endl;
  this->initializeAD(dx.increment());
  Log::trace() << "LinearModelBase<MODEL>::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::stepAD(Increment_ & dx, ModelAuxIncr_ & merr) const {
  Log::trace() << "LinearModelBase<MODEL>::stepAD starting" << std::endl;
  this->stepAD(dx.increment(), merr.modelauxincrement());
  Log::trace() << "LinearModelBase<MODEL>::stepAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::finalizeAD(Increment_ & dx) const {
  Log::trace() << "LinearModelBase<MODEL>::finalizeAD starting" << std::endl;
  this->finalizeAD(dx.increment());
  Log::trace() << "LinearModelBase<MODEL>::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_LINEARMODELBASE_H_
