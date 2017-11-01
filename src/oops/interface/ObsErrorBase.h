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

#include <boost/noncopyable.hpp>
#include <map>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"
#include "util/abor1_cpp.h"
#include "util/Logger.h"
#include "util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Base class for observation error covariance matrices.
/*!
 *  Base class for observation error covariance matrices for a given model.
 *  The interface for the observation error comprises two levels (ObsErrorCovariance
 *  and ObsErrorBase) because we want run time polymorphism.
 *  The ObsErrorCovariance does conversion of arguments to templated ObsVector and
 *  the tracing and timing. The ObsErrorBase does the conversion to model specific
 *  ObsVector.
 */

template<typename MODEL>
class ObsErrorBase : public util::Printable,
                     private boost::noncopyable {
  typedef ObsVector<MODEL>        ObsVector_;

 public:
  ObsErrorBase() {}
  virtual ~ObsErrorBase() {}

/// Linearize and reset for inner loop if needed
  void linearize(const ObsVector_ &);

/// Multiply a Departure by \f$R\f$ and \f$R^{-1}\f$
  ObsVector_ * multiply(const ObsVector_ &) const;
  ObsVector_ * inverseMultiply(const ObsVector_ &) const;

/// Generate random perturbation
  void randomize(ObsVector_ &) const;

/// Get mean error for Jo table
  virtual double getRMSE() const =0;

 private:
  virtual void linearize(const typename MODEL::ObsVector &) =0;
  virtual typename MODEL::ObsVector * multiply(const typename MODEL::ObsVector &) const =0;
  virtual typename MODEL::ObsVector * inverseMultiply(const typename MODEL::ObsVector &) const =0;
  virtual void randomize(typename MODEL::ObsVector &) const =0;
  virtual void print(std::ostream &) const =0;
};

// =============================================================================

/// ObsErrorFactory Factory
template <typename MODEL>
class ObsErrorFactory {
  typedef ObservationSpace<MODEL> ObsSpace_;
 public:
  static ObsErrorBase<MODEL> * create(const ObsSpace_ &, const eckit::Configuration &);
  virtual ~ObsErrorFactory() { getMakers().clear(); }
 protected:
  explicit ObsErrorFactory(const std::string &);
 private:
  virtual ObsErrorBase<MODEL> * make(const ObsSpace_ &, const eckit::Configuration &) =0;
  static std::map < std::string, ObsErrorFactory<MODEL> * > & getMakers() {
    static std::map < std::string, ObsErrorFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class ObsErrorMaker : public ObsErrorFactory<MODEL> {
  typedef ObservationSpace<MODEL> ObsSpace_;
  virtual ObsErrorBase<MODEL> * make(const ObsSpace_ & obs, const eckit::Configuration & conf)
    { return new T(obs.observationspace(), conf); }
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
ObsErrorBase<MODEL>* ObsErrorFactory<MODEL>::create(const ObsSpace_ & obs,
                                                    const eckit::Configuration & conf) {
  Log::trace() << "ObsErrorBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("covariance");
  typename std::map<std::string, ObsErrorFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in observation error factory." << std::endl;
    ABORT("Element does not exist in ObsErrorFactory.");
  }
  ObsErrorBase<MODEL> * ptr = jerr->second->make(obs, conf);
  Log::trace() << "ObsErrorBase<MODEL>::create done" << std::endl;
  return ptr;
}

// =============================================================================

template <typename MODEL>
void ObsErrorBase<MODEL>::linearize(const ObsVector_ & yy) {
  Log::trace() << "ObsErrorBase<MODEL>::linearize starting" << std::endl;
  this->linearize(yy.obsvector());
  Log::trace() << "ObsErrorBase<MODEL>::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsVector<MODEL> * ObsErrorBase<MODEL>::multiply(const ObsVector_ & dy) const {
  Log::trace() << "ObsErrorBase<MODEL>::multiply starting" << std::endl;
  ObsVector_ * dz = new ObsVector_(this->multiply(dy.obsvector()));
  Log::trace() << "ObsErrorBase<MODEL>::multiply done" << std::endl;
  return dz;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsVector<MODEL> * ObsErrorBase<MODEL>::inverseMultiply(const ObsVector_ & dy) const {
  Log::trace() << "ObsErrorBase<MODEL>::inverseMultiply starting" << std::endl;
  ObsVector_ * dz = new ObsVector_(this->inverseMultiply(dy.obsvector()));
  Log::trace() << "ObsErrorBase<MODEL>::inverseMultiply done" << std::endl;
  return dz;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsErrorBase<MODEL>::randomize(ObsVector_ & dy) const {
  Log::trace() << "ObsErrorBase<MODEL>::randomize starting" << std::endl;
  this->randomize(dy.obsvector());
  Log::trace() << "ObsErrorBase<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSERRORBASE_H_
