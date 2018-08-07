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

#include <map>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/generic/VariableChangeBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

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

 public:
  ModelSpaceCovarianceBase(const Geometry_ & resol, const Variables & vars,
                           const eckit::Configuration & conf, const State_ & xb) {
    if (config.has("balance")) {
      eckit::LocalConfiguration balConf("balance", config);
      balop_.reset(VariableChangeFactory<MODEL>::create(resol, vars, balConf, xb));
    }
  }
  virtual ~ModelSpaceCovarianceBase() {}

  const VariableChangeBase_ & getK() const {return *balop_;}
  bool hasK() const {return (balop_ == 0) ? false : true;}

  void linearize(const State_ & fg, const Geometry_ & geom) {
    if (balop_) balop_->doLinearize(fg, geom);
    this->doLinearize(fg, geom);
  }

  void multiply(const Increment_ & dxi, Increment_ & dxo) const {
    if (balop_) {
      Increment_ tmpin = balop_->transformAdjoint(dxi);
      Increment_ tmpout(tmpin);
      this->doMultiply(tmpin, tmpout);
      balop_->transform(tmpout, dxo);
    } else {
      this->doMultiply(dxi, dxo);
    }
  }

  void inverseMultiply(const Increment_ & dxi, Increment_ & dxo) const {
    if (balop_) {
      Increment_ tmp(dxi);
      Increment_ tmpin = balop_->transformInverse(dxi);
      Increment_ tmpout(tmpin);
      this->doInverseMultiply(tmpin, tmpout);
      balop_->transformAdjointInverse(tmpout, dxo);
    } else {
      this->doInverseMultiply(dxi, dxo);
    }
  }

  virtual void randomize(Increment_ &) const = 0;

 private:
  virtual void doLinearize(const State_ &, const Geometry_ &) = 0;
  virtual void doMultiply(const Increment_ &, Increment_ &) const = 0;
  virtual void doInverseMultiply(const Increment_ &, Increment_ &) const = 0;

  boost::scoped_ptr<VariableChangeBase_> balop_;
};

// -----------------------------------------------------------------------------

/// Covariance Factory
template <typename MODEL>
class CovarianceFactory {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;

 public:
  static ModelSpaceCovarianceBase<MODEL> * create(const eckit::Configuration &, const Geometry_ &,
                                                  const Variables &, const State_ &);
  virtual ~CovarianceFactory() { makers_.clear(); }
 protected:
  explicit CovarianceFactory(const std::string &);
 private:
  virtual ModelSpaceCovarianceBase<MODEL> * make(const eckit::Configuration &, const Geometry_ &,
                                                 const Variables &, const State_ &) = 0;
  static std::map < std::string, CovarianceFactory<MODEL> * > makers_;
};

// -----------------------------------------------------------------------------

template<class MODEL, class COVAR>
class CovarMaker : public CovarianceFactory<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;

  virtual ModelSpaceCovarianceBase<MODEL> * make(const eckit::Configuration & conf,
                                                 const Geometry_ & resol,
                                                 const Variables & vars,
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
                                                         const Variables & vars,
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
