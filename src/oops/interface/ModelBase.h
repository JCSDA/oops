/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_MODELBASE_H_
#define OOPS_INTERFACE_MODELBASE_H_

#include <memory>
#include <string>

#include <boost/make_unique.hpp>

#include "oops/base/State.h"
#include "oops/generic/ModelBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/Logger.h"

namespace oops {

namespace interface {

// -----------------------------------------------------------------------------

/// \brief Base class for MODEL-specific implementations of the Model interface.
/// interface::ModelBase overrides oops::ModelBase methods to pass MODEL-specific
/// implementations of State and ModelAuxControl to the MODEL-specific
/// implementation of Model.
///
/// Note: implementations of this interface can opt to extract their settings either from
/// a Configuration object or from a subclass of ModelParametersBase.
///
/// In the former case, they should provide a constructor with the following signature:
///
///    ModelBase(const Geometry_ &, const eckit::Configuration &);
///
/// In the latter case, the implementer should first define a subclass of ModelParametersBase
/// holding the settings of the model in question. The implementation of the ModelBase interface
/// should then typedef `Parameters_` to the name of that subclass and provide a constructor with
/// the following signature:
///
///    ModelBase(const Geometry_ &, const Parameters_ &);
///
template <typename MODEL>
class ModelBase : public oops::ModelBase<MODEL> {
  typedef typename MODEL::ModelAuxControl   ModelAux_;
  typedef typename MODEL::State             State_;

 public:
  static const std::string classname() {return "oops::interface::ModelBase";}

  ModelBase() = default;
  virtual ~ModelBase() = default;

  /// Overrides for oops::ModelBase classes, passing MODEL-specific classes to the
  /// MODEL-specific implementations of Model
  void initialize(oops::State<MODEL> & xx) const final
       { this->initialize(xx.state()); }
  void step(oops::State<MODEL> & xx, const ModelAuxControl<MODEL> & modelaux) const final
       { this->step(xx.state(), modelaux.modelauxcontrol()); }
  void finalize(oops::State<MODEL> & xx) const final
       { this->finalize(xx.state()); }

  /// \brief Forecast initialization, called before every forecast run
  virtual void initialize(State_ &) const = 0;
  /// \brief Forecast "step", called during forecast run; updates state to the next time
  virtual void step(State_ &, const ModelAux_ &) const = 0;
  /// \brief Forecast finalization; called after each forecast run
  virtual void finalize(State_ &) const = 0;
};

// -----------------------------------------------------------------------------

/// \brief A subclass of ModelFactory able to create instances of T (a concrete subclass of
/// interface::ModelBase<MODEL>). Passes MODEL::Geometry to the constructor of T.
template<class MODEL, class T>
class ModelMaker : public ModelFactory<MODEL> {
 private:
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// GenericModelParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<T, GenericModelParameters> Parameters_;

 public:
  typedef oops::Geometry<MODEL>   Geometry_;

  explicit ModelMaker(const std::string & name) : ModelFactory<MODEL>(name) {}

  oops::ModelBase<MODEL> * make(const Geometry_ & geom,
                                const ModelParametersBase & parameters) override {
    Log::trace() << "interface::ModelBase<MODEL>::make starting" << std::endl;
    const auto &stronglyTypedParameters = dynamic_cast<const Parameters_&>(parameters);
    return new T(geom.geometry(),
                 parametersOrConfiguration<HasParameters_<T>::value>(stronglyTypedParameters));
  }

  std::unique_ptr<ModelParametersBase> makeParameters() const override {
    return boost::make_unique<Parameters_>();
  }
};

// -----------------------------------------------------------------------------

}  // namespace interface

}  // namespace oops

#endif  // OOPS_INTERFACE_MODELBASE_H_
