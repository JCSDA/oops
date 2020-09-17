/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_MODELBASE_H_
#define OOPS_BASE_MODELBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/make_unique.hpp>
#include <boost/noncopyable.hpp>

#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
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

/// Base class for encapsulation of the forecast model.
/*!
 * Defines the interfaces for a forecast model.
 */

template <typename MODEL>
class ModelBase : public util::Printable,
                  private boost::noncopyable {
  typedef Geometry<MODEL>            Geometry_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::ModelBase";}

  ModelBase() {}
  virtual ~ModelBase() {}

// Run the model forecast
  void initialize(State_ &) const;
  void step(State_ &, const ModelAux_ &) const;
  void finalize(State_ &) const;

// Information and diagnostics
  virtual const util::Duration & timeResolution() const = 0;
  virtual const oops::Variables & variables() const = 0;

 protected:
// Run the model forecast
  virtual void initialize(typename MODEL::State &) const = 0;
  virtual void step(typename MODEL::State &, const typename MODEL::ModelAuxControl &) const = 0;
  virtual void finalize(typename MODEL::State &) const = 0;

 private:
// Information and diagnostics
  virtual void print(std::ostream &) const = 0;
};

// =============================================================================

template <typename MODEL>
class ModelFactory;

// -----------------------------------------------------------------------------

/// \brief Base class for classes storing model-specific parameters.
class ModelParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(ModelParametersBase, Parameters)
 public:
  /// \brief Model name.
  ///
  /// \note This parameter is marked as optional because it is only required in certain
  /// circumstances (e.g. when model parameters are deserialized into a ModelParametersWrapper
  /// and used by ModelFactory to instantiate a model whose type is determined at runtime), but
  /// not others (e.g. in tests written with a particular model in mind). ModelParametersWrapper
  /// will throw an exception if this parameter is not provided.
  OptionalParameter<std::string> name{"name", this};
};

// -----------------------------------------------------------------------------

/// \brief A subclass of ModelParametersBase storing the values of all options in a
/// single Configuration object.
///
/// This object can be accessed by calling the value() method of the \p config member variable.
///
/// The ConfigurationParameter class does not perform any parameter validation; models using
/// GenericModelParameters should therefore ideally be refactored, replacing this class with a
/// dedicated subclass of ModelParametersBase storing each parameter in a separate
/// (Optional/Required)Parameter object.
class GenericModelParameters : public ModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericModelParameters, ModelParametersBase)
 public:
  ConfigurationParameter config{this};
};

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// ModelParametersBase.
template <typename MODEL>
class ModelParametersWrapper : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelParametersWrapper, Parameters)
 public:
  /// After deserializtion, holds an instance of a subclass of ModelParametersBase controlling
  /// the behavior of a model. The type of the subclass is determined by the value of the "name"
  /// key in the Configuration object from which this object is deserialized.
  RequiredPolymorphicParameter<ModelParametersBase, ModelFactory<MODEL>>
    modelParameters{"name", this};
};

// =============================================================================

/// Model factory
template <typename MODEL>
class ModelFactory {
  typedef Geometry<MODEL>   Geometry_;

 public:
  /// \brief Create and return a new model.
  ///
  /// The model's type is determined by the \c name attribute of \p parameters.
  /// \p parameters must be an instance of the subclass of ModelParametersBase
  /// associated with that model type, otherwise an exception will be thrown.
  static ModelBase<MODEL> * create(const Geometry_ &, const ModelParametersBase & parameters);

  /// \brief Create and return an instance of the subclass of ModelParametersBase
  /// storing parameters of models of the specified type.
  static std::unique_ptr<ModelParametersBase> createParameters(const std::string &name);

  /// \brief Return the names of all models that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~ModelFactory() = default;

 protected:
  /// \brief Register a maker able to create models of type \p name.
  explicit ModelFactory(const std::string & name);

 private:
  virtual ModelBase<MODEL> * make(const Geometry_ &, const ModelParametersBase &) = 0;

  virtual std::unique_ptr<ModelParametersBase> makeParameters() const = 0;

  static std::map < std::string, ModelFactory<MODEL> * > & getMakers() {
    static std::map < std::string, ModelFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

/// \brief A subclass of ModelFactory able to create instances of T (a concrete subclass of
/// ModelBase<MODEL>).
template<class MODEL, class T>
class ModelMaker : public ModelFactory<MODEL> {
 private:
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// GenericModelParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<T, GenericModelParameters> Parameters_;

 public:
  typedef Geometry<MODEL>   Geometry_;

  explicit ModelMaker(const std::string & name) : ModelFactory<MODEL>(name) {}

  ModelBase<MODEL> * make(const Geometry_ & geom, const ModelParametersBase & parameters) override {
    const auto &stronglyTypedParameters = dynamic_cast<const Parameters_&>(parameters);
    return new T(geom.geometry(),
                 parametersOrConfiguration<HasParameters_<T>::value>(stronglyTypedParameters));
  }

  std::unique_ptr<ModelParametersBase> makeParameters() const override {
    return boost::make_unique<Parameters_>();
  }
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelFactory<MODEL>::ModelFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in the model factory."  << std::endl;
    ABORT("Element already registered in ModelFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelBase<MODEL> * ModelFactory<MODEL>::create(const Geometry_ & geom,
                                               const ModelParametersBase & parameters) {
  Log::trace() << "ModelFactory<MODEL>::create starting" << std::endl;
  const std::string &id = parameters.name.value().value();
  typename std::map<std::string, ModelFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in the model factory." << std::endl;
    ABORT("Element does not exist in ModelFactory.");
  }
  ModelBase<MODEL> * ptr = jerr->second->make(geom, parameters);
  Log::trace() << "ModelFactory<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<ModelParametersBase> ModelFactory<MODEL>::createParameters(
    const std::string &name) {
  typename std::map<std::string, ModelFactory<MODEL>*>::iterator it = getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in the model factory");
  }
  return it->second->makeParameters();
}

// =============================================================================

template<typename MODEL>
void ModelBase<MODEL>::initialize(State_ & xx) const {
  Log::trace() << "ModelBase<MODEL>::initialize starting" << std::endl;
  this->initialize(xx.state());
  Log::trace() << "ModelBase<MODEL>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelBase<MODEL>::step(State_ & xx, const ModelAux_ & merr) const {
  Log::trace() << "ModelBase<MODEL>::step starting" << std::endl;
  this->step(xx.state(), merr.modelauxcontrol());
  Log::trace() << "ModelBase<MODEL>::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelBase<MODEL>::finalize(State_ & xx) const {
  Log::trace() << "ModelBase<MODEL>::finalize starting" << std::endl;
  this->finalize(xx.state());
  Log::trace() << "ModelBase<MODEL>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_MODELBASE_H_
