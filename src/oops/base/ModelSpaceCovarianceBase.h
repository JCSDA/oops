/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_MODELSPACECOVARIANCEBASE_H_
#define OOPS_BASE_MODELSPACECOVARIANCEBASE_H_

#include <boost/noncopyable.hpp>
#include <map>
#include <string>

#include "util/Logger.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "eckit/config/Configuration.h"
#include "util/abor1_cpp.h"

namespace util {
  class DateTime;
}

namespace oops {

// -----------------------------------------------------------------------------

// Should this be one with the ErrorCovariance class in the interface directory? YT

/// Abstract base class for model space error covariances.
template <typename MODEL>
class ModelSpaceCovarianceBase {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Increment<MODEL>           Increment_;
  typedef Variables<MODEL>           Variables_;

 public:
  virtual ~ModelSpaceCovarianceBase() {}

  virtual void linearize(const State_ &, const Geometry_ &) =0;
  virtual void multiply(const Increment_ &, Increment_ &) const =0;
  virtual void inverseMultiply(const Increment_ &, Increment_ &) const =0;
  virtual void randomize(Increment_ &) const =0;
};

// -----------------------------------------------------------------------------

/// Covariance Factory
template <typename MODEL>
class CovarianceFactory {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Variables<MODEL>           Variables_;

 public:
  static ModelSpaceCovarianceBase<MODEL> * create(const eckit::Configuration &, const Geometry_ &,
                                                  const Variables_ &, const State_ &);
  virtual ~CovarianceFactory() { makers_.clear(); }
 protected:
  explicit CovarianceFactory(const std::string &);
 private:
  virtual ModelSpaceCovarianceBase<MODEL> * make(const eckit::Configuration &, const Geometry_ &,
                                                 const Variables_ &, const State_ &) =0;
  static std::map < std::string, CovarianceFactory<MODEL> * > makers_;
};

// -----------------------------------------------------------------------------

template<class MODEL, class COVAR>
class CovarMaker : public CovarianceFactory<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Variables<MODEL>           Variables_;

  virtual ModelSpaceCovarianceBase<MODEL> * make(const eckit::Configuration & conf,
                                                 const Geometry_ & resol,
                                                 const Variables_ & vars,
                                                 const State_ & xb) {
    return new COVAR(resol, vars, conf, xb);
  }
 public:
  explicit CovarMaker(const std::string & name) : CovarianceFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
std::map < std::string, CovarianceFactory<MODEL> * > CovarianceFactory<MODEL>::makers_;

// -----------------------------------------------------------------------------

template <typename MODEL>
CovarianceFactory<MODEL>::CovarianceFactory(const std::string & name) {
  if (makers_.find(name) != makers_.end()) {
    Log::error() << name << " already registered in covariance factory." << std::endl;
    ABORT("Element already registered in CovarianceFactory.");
  }
  makers_[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>* CovarianceFactory<MODEL>::create(
                                                         const eckit::Configuration & conf,
                                                         const Geometry_ & resol,
                                                         const Variables_ & vars,
                                                         const State_ & xb) {
  const std::string id = conf.getString("covariance");
  Log::trace() << "ModelSpaceCovarianceBase type = " << id << std::endl;
  typename std::map<std::string, CovarianceFactory<MODEL>*>::iterator
    j = makers_.find(id);
  if (j == makers_.end()) {
    Log::error() << id << " does not exist in covariance factory." << std::endl;
    Log::error() << "Covariance factory contains " << makers_.size() << " elements:" << std::endl;
    for (typename std::map<std::string, CovarianceFactory<MODEL>*>::const_iterator
         jj = makers_.begin(); jj != makers_.end(); ++jj) {
       Log::error() << "A" << jj->first << "B" << std::endl;
    }
    ABORT("Element does not exist in CovarianceFactory.");
  }
  return (*j).second->make(conf, resol, vars, xb);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_MODELSPACECOVARIANCEBASE_H_
