/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_LINEARMODELBASE_H_
#define OOPS_INTERFACE_LINEARMODELBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/interface/State.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

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

/// LinearModelFactory Factory
template <typename MODEL>
class LinearModelFactory {
  typedef Geometry<MODEL>   Geometry_;
 public:
  static LinearModelBase<MODEL> * create(const Geometry_ &, const eckit::Configuration &);
  virtual ~LinearModelFactory() { getMakers().clear(); }
 protected:
  explicit LinearModelFactory(const std::string &);
 private:
  virtual LinearModelBase<MODEL> * make(const Geometry_ &, const eckit::Configuration &) = 0;
  static std::map < std::string, LinearModelFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LinearModelFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LinearModelMaker : public LinearModelFactory<MODEL> {
  typedef Geometry<MODEL>   Geometry_;
  virtual LinearModelBase<MODEL> * make(const Geometry_ & geom, const eckit::Configuration & conf)
    { return new T(geom.geometry(), conf); }
 public:
  explicit LinearModelMaker(const std::string & name) : LinearModelFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearModelFactory<MODEL>::LinearModelFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in the tangent linear model factory." << std::endl;
    ABORT("Element already registered in LinearModelFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearModelBase<MODEL>* LinearModelFactory<MODEL>::create(const Geometry_ & geom,
                                                    const eckit::Configuration & conf) {
  Log::trace() << "LinearModelBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("version");
  typename std::map<std::string, LinearModelFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in the tangent linear model factory." << std::endl;
    ABORT("Element does not exist in LinearModelFactory.");
  }
  LinearModelBase<MODEL> * ptr = jerr->second->make(geom, conf);
  Log::trace() << "LinearModelBase<MODEL>::create done" << std::endl;
  return ptr;
}

// =============================================================================

template<typename MODEL>
void LinearModelBase<MODEL>::setTrajectory(const State_ & xx, State_ & xlr,
                                       const ModelAux_ & maux) {
  Log::trace() << "LinearModelBase<MODEL>::setTrajectory starting" << std::endl;
  util::Timer timer(classname(), "setTrajectory");
  this->setTrajectory(xx.state(), xlr.state(), maux.modelauxcontrol());
  Log::trace() << "LinearModelBase<MODEL>::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::initializeTL(Increment_ & dx) const {
  Log::trace() << "LinearModelBase<MODEL>::initializeTL starting" << std::endl;
  util::Timer timer(classname(), "initializeTL");
  this->initializeTL(dx.increment());
  Log::trace() << "LinearModelBase<MODEL>::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::stepTL(Increment_ & dx, const ModelAuxIncr_ & merr) const {
  Log::trace() << "LinearModelBase<MODEL>::stepTL starting" << std::endl;
  util::Timer timer(classname(), "stepTL");
  this->stepTL(dx.increment(), merr.modelauxincrement());
  Log::trace() << "LinearModelBase<MODEL>::stepTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::finalizeTL(Increment_ & dx) const {
  Log::trace() << "LinearModelBase<MODEL>::finalizeTL starting" << std::endl;
  util::Timer timer(classname(), "finalizeTL");
  this->finalizeTL(dx.increment());
  Log::trace() << "LinearModelBase<MODEL>::finalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::initializeAD(Increment_ & dx) const {
  Log::trace() << "LinearModelBase<MODEL>::initializeAD starting" << std::endl;
  util::Timer timer(classname(), "initializeAD");
  this->initializeAD(dx.increment());
  Log::trace() << "LinearModelBase<MODEL>::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::stepAD(Increment_ & dx, ModelAuxIncr_ & merr) const {
  Log::trace() << "LinearModelBase<MODEL>::stepAD starting" << std::endl;
  util::Timer timer(classname(), "stepAD");
  this->stepAD(dx.increment(), merr.modelauxincrement());
  Log::trace() << "LinearModelBase<MODEL>::stepAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModelBase<MODEL>::finalizeAD(Increment_ & dx) const {
  Log::trace() << "LinearModelBase<MODEL>::finalizeAD starting" << std::endl;
  util::Timer timer(classname(), "finalizeAD");
  this->finalizeAD(dx.increment());
  Log::trace() << "LinearModelBase<MODEL>::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEARMODELBASE_H_
