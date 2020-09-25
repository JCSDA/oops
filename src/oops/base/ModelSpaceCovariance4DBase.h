/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_MODELSPACECOVARIANCE4DBASE_H_
#define OOPS_BASE_MODELSPACECOVARIANCE4DBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
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
class ModelSpaceCovariance4DBase {
  typedef Geometry<MODEL>                                       Geometry_;
  typedef State4D<MODEL>                                        State4D_;
  typedef Increment<MODEL>                                      Increment_;
  typedef Increment4D<MODEL>                                    Increment4D_;
  typedef LinearVariableChangeBase<MODEL>                       LinearVariableChangeBase_;
  typedef typename boost::ptr_vector<LinearVariableChangeBase_> ChvarVec_;
  typedef typename ChvarVec_::iterator iter_;
  typedef typename ChvarVec_::const_iterator icst_;
  typedef typename ChvarVec_::const_reverse_iterator ircst_;

 public:
  ModelSpaceCovariance4DBase(const State4D_ &, const State4D_ &,
                             const Geometry_ &, const eckit::Configuration &);
  virtual ~ModelSpaceCovariance4DBase() {}

//  const LinearVariableChangeBase_ & getK(const unsigned & ii) const {return *chvars_[ii];}
//  bool hasK() const { return (chvars_.size() == 0) ? false : true; }

  void randomize(Increment4D_ &) const;
  void multiply(const Increment4D_ &, Increment4D_ &) const;
  void inverseMultiply(const Increment4D_ &, Increment4D_ &) const;

 private:
  virtual void doRandomize(Increment4D_ &) const = 0;
  virtual void doMultiply(const Increment4D_ &, Increment4D_ &) const = 0;
  virtual void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const = 0;

  ChvarVec_ chvars_;
};

// -----------------------------------------------------------------------------

/// Covariance Factory
template <typename MODEL>
class Covariance4DFactory {
  typedef Geometry<MODEL>            Geometry_;
  typedef State4D<MODEL>             State4D_;

 public:
  static ModelSpaceCovariance4DBase<MODEL> * create(const eckit::Configuration &,
                                                    const Geometry_ &, const Variables &,
                                                    const State4D_ &, const State4D_ &);
  virtual ~Covariance4DFactory() = default;
 protected:
  explicit Covariance4DFactory(const std::string &);
 private:
  virtual ModelSpaceCovariance4DBase<MODEL> * make(const eckit::Configuration &,
                                                   const Geometry_ &, const Variables &,
                                                   const State4D_ &, const State4D_ &) = 0;
  static std::map < std::string, Covariance4DFactory<MODEL> * > & getMakers() {
    static std::map < std::string, Covariance4DFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class COVAR>
class Covar4DMaker : public Covariance4DFactory<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef State4D<MODEL>             State4D_;

  virtual ModelSpaceCovariance4DBase<MODEL> * make(const eckit::Configuration & conf,
                                                   const Geometry_ & resol, const Variables & vars,
                                                   const State4D_ & xb, const State4D_ & fg) {
    return new COVAR(resol, vars, conf, xb, fg);
  }
 public:
  explicit Covar4DMaker(const std::string & name) : Covariance4DFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
Covariance4DFactory<MODEL>::Covariance4DFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in covariance factory." << std::endl;
    ABORT("Element already registered in Covariance4DFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelSpaceCovariance4DBase<MODEL>* Covariance4DFactory<MODEL>::create(
                                                         const eckit::Configuration & conf,
                                                         const Geometry_ & resol,
                                                         const Variables & vars,
                                                         const State4D_ & xb, const State4D_ & fg) {
  const std::string id = conf.getString("covariance model");
  Log::trace() << "ModelSpaceCovariance4DBase type = " << id << std::endl;
  typename std::map<std::string, Covariance4DFactory<MODEL>*>::iterator jcov = getMakers().find(id);
  if (jcov == getMakers().end()) {
    Log::error() << id << " does not exist in Covariance4DFactory." << std::endl;
    Log::error() << "Covariance4DFactory has " << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, Covariance4DFactory<MODEL>*>::const_iterator
         jj = getMakers().begin(); jj != getMakers().end(); ++jj) {
       Log::error() << "A " << jj->first << " B" << std::endl;
    }
    ABORT("Element does not exist in Covariance4DFactory.");
  }
  Variables vars_in(vars);
  Variables vars_out;
  if (conf.has("variable changes")) {
    std::vector<eckit::LocalConfiguration> chvarconfs;
    conf.get("variable changes", chvarconfs);
    for (const auto & config : boost::adaptors::reverse(chvarconfs)) {
      vars_out = Variables(config, "output variables");
      if (!(vars_in == vars_out)) {
        Log::error() << "Input variables:  " << vars_in << std::endl;
        Log::error() << "Output variables: " << vars_out << std::endl;
        ABORT("Sequence of variable changes is not consistent");
      }
      vars_in = Variables(config, "input variables");
    }
  }
  return (*jcov).second->make(conf, resol, vars_in, xb, fg);
}

// =============================================================================

template <typename MODEL>
ModelSpaceCovariance4DBase<MODEL>::ModelSpaceCovariance4DBase(const State4D_ & bg,
                                                              const State4D_ & fg,
                                                              const Geometry_ & resol,
                                                              const eckit::Configuration & conf) {
  ASSERT(bg.size() == fg.size());
  if (conf.has("variable changes")) {
    std::vector<eckit::LocalConfiguration> chvarconfs;
    conf.get("variable changes", chvarconfs);
    for (const auto & config : chvarconfs) {
      chvars_.push_back(LinearVariableChangeFactory<MODEL>::create(bg[0], fg[0], resol, config));
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovariance4DBase<MODEL>::randomize(Increment4D_ & dx) const {
  // TODO(notguillaume): Generalize to non-square change of variable
  if (chvars_.size()) {
    this->doRandomize(dx);   // dx = C^1/2 dx
    // K_N K_N-1 ... K_1
    std::unique_ptr<Increment4D_> dxchvarin4D(new Increment4D_(dx, true));
    for (int isub = dx.first(); isub <= dx.last(); ++isub) {
      std::unique_ptr<Increment_> dxchvarin(new Increment_((*dxchvarin4D)[isub], true));
      for (icst_ it = chvars_.begin(); it != chvars_.end(); ++it) {
        Increment_ dxchvarout = it->multiply(*dxchvarin);
        dxchvarin.reset(new Increment_(dxchvarout));
      }
      (*dxchvarin4D)[isub] = *dxchvarin;
    }
    dx = *dxchvarin4D;
  } else {
    this->doRandomize(dx);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovariance4DBase<MODEL>::multiply(const Increment4D_ & dxi,
                                               Increment4D_ & dxo) const {
  ASSERT(dxi.last() == dxo.last());
  if (chvars_.size()) {
    // K_1^T K_2^T .. K_N^T
    std::unique_ptr<Increment4D_> dxchvarin4D(new Increment4D_(dxi, true));
    for (int isub = dxi.first(); isub <= dxi.last(); ++isub) {
      std::unique_ptr<Increment_> dxchvarin(new Increment_((*dxchvarin4D)[isub], true));
      for (ircst_ it = chvars_.rbegin(); it != chvars_.rend(); ++it) {
        Increment_ dxchvarout = it->multiplyAD(*dxchvarin);
        dxchvarin.reset(new Increment_(dxchvarout));
      }
      (*dxchvarin4D)[isub] = *dxchvarin;
    }
    Increment4D_ dxchvarout4D(*dxchvarin4D, false);

    this->doMultiply(*dxchvarin4D, dxchvarout4D);

    // K_N K_N-1 ... K_1
    dxchvarin4D.reset(new Increment4D_(dxchvarout4D));
    for (int isub = dxi.first(); isub <= dxi.last(); ++isub) {
      std::unique_ptr<Increment_> dxchvarin(new Increment_((*dxchvarin4D)[isub], true));
      for (icst_ it = chvars_.begin(); it != chvars_.end(); ++it) {
        Increment_ dxchvarout = it->multiply(*dxchvarin);
        dxchvarin.reset(new Increment_(dxchvarout));
      }
      (*dxchvarin4D)[isub] = *dxchvarin;
    }
    dxo = *dxchvarin4D;
  } else {
    this->doMultiply(dxi, dxo);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovariance4DBase<MODEL>::inverseMultiply(const Increment4D_ & dxi,
                                                        Increment4D_ & dxo) const {
  ASSERT(dxi.last() == dxo.last());
  if (chvars_.size()) {
    // K_1^{-1} K_2^{-1} .. K_N^{-1}
    std::unique_ptr<Increment4D_> dxchvarin4D(new Increment4D_(dxi, true));
    for (int isub = dxi.first(); isub <= dxi.last(); ++isub) {
      std::unique_ptr<Increment_> dxchvarin(new Increment_((*dxchvarin4D)[isub], true));
      for (ircst_ it = chvars_.rbegin(); it != chvars_.rend(); ++it) {
        Increment_ dxchvarout = it->multiplyInverse(*dxchvarin);
        dxchvarin.reset(new Increment_(dxchvarout));
      }
      (*dxchvarin4D)[isub] = *dxchvarin;
    }
    Increment4D_ dxchvarout4D(*dxchvarin4D, false);

    this->doInverseMultiply(*dxchvarin4D, dxchvarout4D);

    // K_N^T^{-1} K_N-1^T^{-1} ... K_1^T^{-1}
    dxchvarin4D.reset(new Increment4D_(dxchvarout4D));
    for (int isub = dxi.first(); isub <= dxi.last(); ++isub) {
      std::unique_ptr<Increment_> dxchvarin(new Increment_((*dxchvarin4D)[isub], true));
      for (icst_ it = chvars_.begin(); it != chvars_.end(); ++it) {
        Increment_ dxchvarout = it->multiplyInverseAD(*dxchvarin);
        dxchvarin.reset(new Increment_(dxchvarout));
      }
      (*dxchvarin4D)[isub] = *dxchvarin;
    }
    dxo = *dxchvarin4D;
  } else {
    this->doInverseMultiply(dxi, dxo);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_MODELSPACECOVARIANCE4DBASE_H_
