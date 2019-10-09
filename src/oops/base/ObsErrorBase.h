/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERRORBASE_H_
#define OOPS_BASE_OBSERRORBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>
#include "eckit/config/Configuration.h"

#include "oops/interface/ObsSpace.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Base class for observation error covariance matrices.

template<typename MODEL>
class ObsErrorBase : public util::Printable,
                     private boost::noncopyable {
  typedef ObsVector<MODEL>        ObsVector_;
  typedef ObsSpace<MODEL>         ObsSpace_;

 public:
  ObsErrorBase() {}
  virtual ~ObsErrorBase() {}

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
  typedef ObsSpace<MODEL> ObsSpace_;
 public:
  static ObsErrorBase<MODEL> * create(const eckit::Configuration &, const ObsSpace_ &);
  virtual ~ObsErrorFactory() { getMakers().clear(); }
 protected:
  explicit ObsErrorFactory(const std::string &);
 private:
  virtual ObsErrorBase<MODEL> * make(const eckit::Configuration &, const ObsSpace_ &) = 0;
  static std::map < std::string, ObsErrorFactory<MODEL> * > & getMakers() {
    static std::map < std::string, ObsErrorFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class ObsErrorMaker : public ObsErrorFactory<MODEL> {
  typedef ObsSpace<MODEL> ObsSpace_;
  virtual ObsErrorBase<MODEL> * make(const eckit::Configuration & conf,
                                     const ObsSpace_ & obs)
    { return new T(conf, obs); }
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
                                                    const ObsSpace_ & obs) {
  Log::trace() << "ObsErrorBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("covariance");
  typename std::map<std::string, ObsErrorFactory<MODEL>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    Log::error() << id << " does not exist in observation error factory." << std::endl;
    ABORT("Element does not exist in ObsErrorFactory.");
  }
  ObsErrorBase<MODEL> * ptr = jerr->second->make(conf, obs);
  Log::trace() << "ObsErrorBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERRORBASE_H_
