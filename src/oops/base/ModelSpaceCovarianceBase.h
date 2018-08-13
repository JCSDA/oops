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
#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
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
  typedef LinearVariableChangeBase<MODEL>  LinearVariableChangeBase_;
  typedef typename boost::ptr_vector<LinearVariableChangeBase_> ChvarVec_;
  typedef typename ChvarVec_::iterator iter_;
  typedef typename ChvarVec_::const_iterator icst_;
  typedef typename ChvarVec_::const_reverse_iterator ircst_;

 public:
  explicit ModelSpaceCovarianceBase(const eckit::Configuration & conf);
  virtual ~ModelSpaceCovarianceBase() {}

  const LinearVariableChangeBase_ & getK(const unsigned & ii) const {return *chvars_[ii];}
  bool hasK() const { return (chvars_.size() == 0) ? false : true; }

  void linearize(const State_ &, const Geometry_ &);
  void multiply(const Increment_ &, Increment_ &) const;
  void inverseMultiply(const Increment_ &, Increment_ &) const;

  virtual void randomize(Increment_ &) const = 0;

 private:
  virtual void doLinearize(const State_ &, const Geometry_ &) = 0;
  virtual void doMultiply(const Increment_ &, Increment_ &) const = 0;
  virtual void doInverseMultiply(const Increment_ &, Increment_ &) const = 0;

  ChvarVec_ chvars_;
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

// =============================================================================

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>::ModelSpaceCovarianceBase(const eckit::Configuration & conf) {
  if (conf.has("variable_changes")) {
    std::vector<eckit::LocalConfiguration> chvarconfs;
    conf.get("variable_changes", chvarconfs);
    for (const auto & config : chvarconfs)
      chvars_.push_back(LinearVariableChangeFactory<MODEL>::create(config));
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovarianceBase<MODEL>::linearize(const State_ & fg,
                                                const Geometry_ & geom) {
  if (chvars_.size()) {
    for (iter_ it = chvars_.begin(); it != chvars_.end(); ++it)
      it->linearize(fg, geom);
  }

  this->doLinearize(fg, geom);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovarianceBase<MODEL>::multiply(const Increment_ & dxi,
                                               Increment_ & dxo) const {
  if (chvars_.size()) {
    // K_1^T K_2^T .. K_N^T
    boost::scoped_ptr<Increment_> dxchvarin(new Increment_(dxi));
    for (ircst_ it = chvars_.rbegin(); it != chvars_.rend(); ++it) {
      Increment_ dxchvarout = it->multiplyAD(*dxchvarin);
      dxchvarin.reset(new Increment_(dxchvarout));
    }
    Increment_ dxchvarout(*dxchvarin, false);

    this->doMultiply(*dxchvarin, dxchvarout);

    // K_N K_N-1 ... K_1
    dxchvarin.reset(new Increment_(dxchvarout));
    for (icst_ it = chvars_.begin(); it != chvars_.end(); ++it) {
      Increment_ dxchvarout = it->multiply(*dxchvarin);
      dxchvarin.reset(new Increment_(dxchvarout));
    }
    dxo = *dxchvarin;
  } else {
    this->doMultiply(dxi, dxo);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovarianceBase<MODEL>::inverseMultiply(const Increment_ & dxi,
                                                      Increment_ & dxo) const {
  if (chvars_.size()) {
    // K_1^{-1} K_2^{-1} .. K_N^{-1}
    boost::scoped_ptr<Increment_> dxchvarin(new Increment_(dxi));
    for (ircst_ it = chvars_.rbegin(); it != chvars_.rend(); ++it) {
      Increment_ dxchvarout = it->multiplyInverse(*dxchvarin);
      dxchvarin.reset(new Increment_(dxchvarout));
    }
    Increment_ dxchvarout(*dxchvarin, false);

    this->doInverseMultiply(*dxchvarin, dxchvarout);

    // K_N^T^{-1} K_N-1^T^{-1} ... K_1^T^{-1}
    dxchvarin.reset(new Increment_(dxchvarout));
    for (icst_ it = chvars_.begin(); it != chvars_.end(); ++it) {
      Increment_ dxchvarout = it->multiplyInverseAD(*dxchvarin);
      dxchvarin.reset(new Increment_(dxchvarout));
    }
    dxo = *dxchvarin;
  } else {
    this->doInverseMultiply(dxi, dxo);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_MODELSPACECOVARIANCEBASE_H_
