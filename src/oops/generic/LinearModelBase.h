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

#include <boost/noncopyable.hpp>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Base class for generic implementations of the linearized forecasting models.
/// Use this class as a base class for generic implementations,
/// and interface::LinearModelBase as a base class for MODEL-specific implementations.
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

/// LinearModel factory
template <typename MODEL>
class LinearModelFactory {
  typedef Geometry<MODEL>   Geometry_;

 public:
  /// \brief Create and return a new linear model.
  static LinearModelBase<MODEL> * create(const Geometry_ &, const eckit::Configuration &);

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
  virtual LinearModelBase<MODEL> * make(const Geometry_ &, const eckit::Configuration &) = 0;

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
 public:
  typedef Geometry<MODEL>   Geometry_;

  explicit LinearModelMaker(const std::string & name) : LinearModelFactory<MODEL>(name) {}

  LinearModelBase<MODEL> * make(const Geometry_ & geom,
                                const eckit::Configuration & config) override {
    Log::trace() << "LinearModelBase<MODEL>::make starting" << std::endl;
    return new T(geom, config);
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
                                                           const eckit::Configuration & config) {
  Log::trace() << "LinearModelFactory<MODEL>::create starting" << std::endl;
  const std::string id = config.getString("name");
  typename std::map<std::string, LinearModelFactory<MODEL>*>::iterator jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in the linear model factory");
  }
  LinearModelBase<MODEL> * ptr = jerr->second->make(geom, config);
  Log::trace() << "LinearModelFactory<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_LINEARMODELBASE_H_
