/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_MODELBASE_H_
#define OOPS_GENERIC_MODELBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/make_unique.hpp>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Base class for generic implementations of the forecasting models.
/// Use this class as a base class for generic implementations,
/// and interface::ModelBase as a base class for MODEL-specific implementations.

template <typename MODEL>
class ModelBase : public util::Printable,
                  private boost::noncopyable {
  typedef ModelAuxControl<MODEL>   ModelAux_;
  typedef State<MODEL>             State_;

 public:
  static const std::string classname() {return "oops::ModelBase";}

  ModelBase() = default;
  virtual ~ModelBase() = default;

  /// \brief Forecast initialization, called before every forecast run
  virtual void initialize(State_ &) const = 0;
  /// \brief Forecast "step", called during forecast run; updates state to the next time
  virtual void step(State_ &, const ModelAux_ &) const = 0;
  /// \brief Forecast finalization; called after each forecast run
  virtual void finalize(State_ &) const = 0;

  /// \brief Time step for running Model's forecast in oops (frequency with which the
  /// State will be updated)
  virtual const util::Duration & timeResolution() const = 0;

 private:
  /// \brief Print; used for logging
  virtual void print(std::ostream &) const = 0;
};

// =============================================================================

/// Model factory
template <typename MODEL>
class ModelFactory {
  typedef Geometry<MODEL>   Geometry_;

 public:
  /// \brief Create and return a new model.
  static ModelBase<MODEL> * create(const Geometry_ &, const eckit::Configuration &);

  /// \brief Return the names of all models that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~ModelFactory() = default;

 protected:
  /// \brief Register a maker able to create models of type \p name.
  explicit ModelFactory(const std::string & name);

 private:
  virtual ModelBase<MODEL> * make(const Geometry_ &, const eckit::Configuration &) = 0;

  static std::map < std::string, ModelFactory<MODEL> * > & getMakers() {
    static std::map < std::string, ModelFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

/// \brief A subclass of ModelFactory able to create instances of T (a concrete subclass of
/// ModelBase<MODEL>). Passes Geometry<MODEL> to the constructor of T.
template<class MODEL, class T>
class ModelMaker : public ModelFactory<MODEL> {
 public:
  typedef Geometry<MODEL>   Geometry_;

  explicit ModelMaker(const std::string & name) : ModelFactory<MODEL>(name) {}

  ModelBase<MODEL> * make(const Geometry_ & geom,
                          const eckit::Configuration & config) override {
    Log::trace() << "ModelBase<MODEL>::make starting" << std::endl;
    return new T(geom, config);
  }
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelFactory<MODEL>::ModelFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in the model factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelBase<MODEL> * ModelFactory<MODEL>::create(const Geometry_ & geom,
                                               const eckit::Configuration & config) {
  Log::trace() << "ModelFactory<MODEL>::create starting" << std::endl;
  const std::string id = config.getString("name");
  typename std::map<std::string, ModelFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in the model factory");
  }
  ModelBase<MODEL> * ptr = jerr->second->make(geom, config);
  Log::trace() << "ModelFactory<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_MODELBASE_H_
