/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_MODELBASE_H_
#define OOPS_BASE_MODELBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
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

/// ModelFactory Factory
template <typename MODEL>
class ModelFactory {
  typedef Geometry<MODEL>   Geometry_;
 public:
  static ModelBase<MODEL> * create(const Geometry_ &, const eckit::Configuration &);
  virtual ~ModelFactory() = default;
 protected:
  explicit ModelFactory(const std::string &);
 private:
  virtual ModelBase<MODEL> * make(const Geometry_ &, const eckit::Configuration &) = 0;
  static std::map < std::string, ModelFactory<MODEL> * > & getMakers() {
    static std::map < std::string, ModelFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class ModelMaker : public ModelFactory<MODEL> {
  typedef Geometry<MODEL>   Geometry_;
  virtual ModelBase<MODEL> * make(const Geometry_ & geom, const eckit::Configuration & conf)
    { return new T(geom.geometry(), conf); }
 public:
  explicit ModelMaker(const std::string & name) : ModelFactory<MODEL>(name) {}
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
                                         const eckit::Configuration & conf) {
  Log::trace() << "ModelBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("name");
  typename std::map<std::string, ModelFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in the model factory." << std::endl;
    ABORT("Element does not exist in ModelFactory.");
  }
  ModelBase<MODEL> * ptr = jerr->second->make(geom, conf);
  Log::trace() << "ModelBase<MODEL>::create done" << std::endl;
  return ptr;
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
