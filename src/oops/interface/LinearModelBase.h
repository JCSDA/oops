/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_LINEARMODELBASE_H_
#define OOPS_INTERFACE_LINEARMODELBASE_H_

#include <memory>
#include <string>

#include <boost/make_unique.hpp>

#include "oops/generic/LinearModelBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"

namespace oops {

namespace interface {

// -----------------------------------------------------------------------------

/// \brief Base class for MODEL-specific implementations of the LinearModel interface.
/// interface::LinearModelBase overrides oops::LinearModelBase methods to pass MODEL-specific
/// implementations of State, Increment and ModelAuxIncrement to the MODEL-specific implementation
/// of LinearModel.
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
class LinearModelBase : public oops::LinearModelBase<MODEL> {
  typedef typename MODEL::Increment         Increment_;
  typedef typename MODEL::ModelAuxControl   ModelAuxCtl_;
  typedef typename MODEL::ModelAuxIncrement ModelAuxInc_;
  typedef typename MODEL::State             State_;

 public:
  static const std::string classname() {return "oops::interface::LinearModelBase";}

  LinearModelBase() = default;
  virtual ~LinearModelBase() = default;

  /// Overrides for oops::LinearModelBase classes, passing MODEL-specific classes to the
  /// MODEL-specific implementations of LinearModel
  void initializeTL(oops::Increment<MODEL> & dx) const final
       { this->initializeTL(dx.increment()); }
  void stepTL(oops::Increment<MODEL> & dx, const ModelAuxIncrement<MODEL> & modelaux) const final
       { this->stepTL(dx.increment(), modelaux.modelauxincrement()); }
  void finalizeTL(oops::Increment<MODEL> & dx) const final
       { this->finalizeTL(dx.increment()); }

  void initializeAD(oops::Increment<MODEL> & dx) const final
      { this->initializeAD(dx.increment()); }
  void stepAD(oops::Increment<MODEL> & dx, ModelAuxIncrement<MODEL> & modelaux) const final
      { this->stepAD(dx.increment(), modelaux.modelauxincrement()); }
  void finalizeAD(oops::Increment<MODEL> & dx) const final
      { this->finalizeAD(dx.increment()); }

  void setTrajectory(const oops::State<MODEL> & xx, oops::State<MODEL> & xxtraj,
                     const ModelAuxControl<MODEL> & modelaux) final
     { this->setTrajectory(xx.state(), xxtraj.state(), modelaux.modelauxcontrol()); }

  /// \brief Tangent linear forecast initialization, called before every run
  virtual void initializeTL(Increment_ &) const = 0;
  /// \brief Tangent linear forecast "step", called during run; updates Increment to the next time
  virtual void stepTL(Increment_ &, const ModelAuxInc_ &) const = 0;
  /// \brief Tangent linear forecast finalization; called after each run
  virtual void finalizeTL(Increment_ &) const = 0;

  /// \brief Adjoint forecast initialization, called before every run
  virtual void initializeAD(Increment_ &) const = 0;
  /// \brief Adjoint forecast "step", called during run; updates increment to the previous time
  virtual void stepAD(Increment_ &, ModelAuxInc_ &) const = 0;
  /// \brief Adjoint forecast finalization; called after each run
  virtual void finalizeAD(Increment_ &) const = 0;

  /// \brief Set the trajectory for the linear model, called after each step of the forecast.
  /// The incoming State is output from the nonlinear forecast. The adjustable State is
  /// interpolated to the resolution of the linear model.
  virtual void setTrajectory(const State_ &, State_ &, const ModelAuxCtl_ &) = 0;

  /// \brief Print, used in logging
  void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

/// \brief A subclass of LinearModelFactory able to create instances of T (a concrete subclass of
/// interface::LinearModelBase<MODEL>). Passes MODEL::Geometry to the constructor of T.
template<class MODEL, class T>
class LinearModelMaker : public LinearModelFactory<MODEL> {
 private:
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// GenericLinearModelParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<T, GenericLinearModelParameters> Parameters_;

 public:
  typedef oops::Geometry<MODEL>   Geometry_;

  explicit LinearModelMaker(const std::string & name) : LinearModelFactory<MODEL>(name) {}

  oops::LinearModelBase<MODEL> * make(const Geometry_ & geom,
                                      const LinearModelParametersBase & parameters) override {
    Log::trace() << "interface::LinearModelBase<MODEL>::make starting" << std::endl;
    const auto &stronglyTypedParameters = dynamic_cast<const Parameters_&>(parameters);
    return new T(geom.geometry(),
                 parametersOrConfiguration<HasParameters_<T>::value>(stronglyTypedParameters));
  }

  std::unique_ptr<LinearModelParametersBase> makeParameters() const override {
    return boost::make_unique<Parameters_>();
  }
};

// -----------------------------------------------------------------------------

}  // namespace interface

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEARMODELBASE_H_
