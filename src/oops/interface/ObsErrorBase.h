/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_OBSERRORBASE_H_
#define OOPS_INTERFACE_OBSERRORBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>
#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Base class for observation error covariance matrices.
/*!
 *  Base class for observation error covariance matrices.
 *  The interface for the observation error comprises two levels (ObsErrorCovariance
 *  and ObsErrorBase) because we want run time polymorphism.
 *  The ObsErrorCovariance does conversion of arguments to templated ObsVector and
 *  the tracing and timing. The ObsErrorBase does the conversion to model specific
 *  ObsVector.
 */

template<typename MODEL>
class ObsErrorBase : public util::Printable,
                     private boost::noncopyable {
  typedef typename MODEL::ObsVector        ObsVector_;
  typedef typename MODEL::ObsSpace         ObsSpace_;

 public:
  ObsErrorBase() {}
  virtual ~ObsErrorBase() {}

/// Update after obs filters
  virtual void update() = 0;

/// Multiply a Departure by \f$R\f$ and \f$R^{-1}\f$
  virtual void multiply(ObsVector_ &) const = 0;
  virtual void inverseMultiply(ObsVector_ &) const = 0;

/// Generate random perturbation
  virtual void randomize(ObsVector_ &) const = 0;

/// Get mean error for Jo table
  virtual double getRMSE() const = 0;
};

// =============================================================================

/// ObsErrorFactory Factory
template <typename MODEL>
class ObsErrorFactory {
  typedef ObservationSpace<MODEL> ObsSpace_;
 public:
  static ObsErrorBase<MODEL> * create(const eckit::Configuration &, const ObsSpace_ &,
                                      const Variables &);
  virtual ~ObsErrorFactory() { getMakers().clear(); }
 protected:
  explicit ObsErrorFactory(const std::string &);
 private:
  virtual ObsErrorBase<MODEL> * make(const eckit::Configuration &, const ObsSpace_ &,
                                     const Variables &) = 0;
  static std::map < std::string, ObsErrorFactory<MODEL> * > & getMakers() {
    static std::map < std::string, ObsErrorFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class ObsErrorMaker : public ObsErrorFactory<MODEL> {
  typedef ObservationSpace<MODEL> ObsSpace_;
  virtual ObsErrorBase<MODEL> * make(const eckit::Configuration & conf, const ObsSpace_ & obs,
                                     const Variables & vars)
    { return new T(conf, obs.observationspace(), vars); }
 public:
  explicit ObsErrorMaker(const std::string & name) : ObsErrorFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
ObsErrorFactory<MODEL>::ObsErrorFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in observation error factory." << std::endl;
    ABORT("Element already registered in ObsErrorFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsErrorBase<MODEL>* ObsErrorFactory<MODEL>::create(const eckit::Configuration & conf,
                                                    const ObsSpace_ & obs, const Variables & vars) {
  Log::trace() << "ObsErrorBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("covariance");
  typename std::map<std::string, ObsErrorFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in observation error factory." << std::endl;
    ABORT("Element does not exist in ObsErrorFactory.");
  }
  ObsErrorBase<MODEL> * ptr = jerr->second->make(conf, obs, vars);
  Log::trace() << "ObsErrorBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSERRORBASE_H_
