/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_MODELSPACECOVARIANCEBASE_H_
#define OOPS_BASE_MODELSPACECOVARIANCEBASE_H_

#include <algorithm>
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include <boost/make_unique.hpp>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/base/Geometry.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"
#include "oops/util/Timer.h"

namespace util {
  class DateTime;
}

namespace oops {

// -----------------------------------------------------------------------------

// Should this be one with the ErrorCovariance class in the interface directory? YT

/// Abstract base class for model space error covariances.
template <typename MODEL>
class ModelSpaceCovarianceBase {
  typedef Geometry<MODEL>                                       Geometry_;
  typedef GeometryIterator<MODEL>                               GeometryIterator_;
  typedef State4D<MODEL>                                        State4D_;
  typedef Increment<MODEL>                                      Increment_;
  typedef Increment4D<MODEL>                                    Increment4D_;
  typedef LinearVariableChange<MODEL>                           LinearVariableChange_;

 public:
  ModelSpaceCovarianceBase(const Geometry_ &, const eckit::Configuration &,
                           const State4D_ &, const State4D_ &);
  virtual ~ModelSpaceCovarianceBase() {}

  void randomize(Increment4D_ &) const;
  void multiply(const Increment4D_ &, Increment4D_ &) const;
  void inverseMultiply(const Increment4D_ &, Increment4D_ &) const;
  void getVariance(Increment4D_ &) const;

  const std::string covarianceModel() const {return covarianceModel_;}
  size_t randomizationSize() const {return randomizationSize_;}

 protected:
  std::unique_ptr<Variables> BVars_;

 private:
  virtual void doRandomize(Increment4D_ &) const = 0;
  virtual void doMultiply(const Increment4D_ &, Increment4D_ &) const = 0;
  virtual void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const = 0;

  std::string covarianceModel_;
  size_t randomizationSize_;
  bool fullInverse_ = false;
  int fullInverseIterations_;
  double fullInverseAccuracy_;
  std::vector<std::unique_ptr<LinearVariableChange_>> linVarChg_;
  std::vector<std::unique_ptr<LinearVariableChange_>> invLinVarChg_;
  std::unique_ptr<Variables> anaVars_;
  std::string timername_;
};

// =============================================================================

/// Covariance Factory
template <typename MODEL>
class CovarianceFactory {
  typedef Geometry<MODEL> Geometry_;
  typedef State4D<MODEL>  State4D_;

 public:
  static ModelSpaceCovarianceBase<MODEL> * create(const Geometry_ &, const Variables &,
                                                  const eckit::Configuration &,
                                                  const State4D_ &, const State4D_ &);

  virtual ~CovarianceFactory() = default;

 protected:
  explicit CovarianceFactory(const std::string &);

 private:
  virtual ModelSpaceCovarianceBase<MODEL> * make(const Geometry_ &, const Variables &,
                                                 const eckit::Configuration &,
                                                 const State4D_ &, const State4D_ &) = 0;

  static std::map < std::string, CovarianceFactory<MODEL> * > & getMakers() {
    static std::map < std::string, CovarianceFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class COVAR>
class CovarMaker : public CovarianceFactory<MODEL> {
  typedef Geometry<MODEL> Geometry_;
  typedef State4D<MODEL>  State4D_;

  virtual ModelSpaceCovarianceBase<MODEL> * make(const Geometry_ & resol, const Variables & vars,
                                                 const eckit::Configuration & conf,
                                                 const State4D_ & xb, const State4D_ & fg) {
    return new COVAR(resol, vars, conf, xb, fg);
  }
 public:
  explicit CovarMaker(const std::string & name) : CovarianceFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
CovarianceFactory<MODEL>::CovarianceFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in covariance factory." << std::endl;
    ABORT("Element already registered in CovarianceFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>* CovarianceFactory<MODEL>::create(const Geometry_ & resol,
                                                         const Variables & vars,
                                                         const eckit::Configuration & conf,
                                                         const State4D_ & xb, const State4D_ & fg)
{
  const std::string id = conf.getString("covariance model");
  Log::trace() << "ModelSpaceCovarianceBase type = " << id << std::endl;
  typename std::map<std::string, CovarianceFactory<MODEL>*>::iterator jcov = getMakers().find(id);
  if (jcov == getMakers().end()) {
    Log::error() << id << " does not exist in CovarianceFactory." << std::endl;
    Log::error() << "CovarianceFactory has " << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, CovarianceFactory<MODEL>*>::const_iterator
         jj = getMakers().begin(); jj != getMakers().end(); ++jj) {
       Log::error() << "A " << jj->first << " B" << std::endl;
    }
    ABORT("Element does not exist in CovarianceFactory.");
  }
  Variables vars_in(vars);
  if (conf.has("linear variable change")) {
    const eckit::LocalConfiguration config(conf, "linear variable change");
    if (config.has("input variables")) {
      vars_in = Variables(config, "input variables");
    }
  }
  return (*jcov).second->make(resol, vars_in, conf, xb, fg);
}

// =============================================================================

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>::ModelSpaceCovarianceBase(const Geometry_ & resol,
                                                          const eckit::Configuration & conf,
                                                          const State4D_ & xb, const State4D_ & fg)
{
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::ModelSpaceCovarianceBase starting" << std::endl;
  covarianceModel_ = conf.getString("covariance model", "none");
  timername_ = "oops::Covariance::" + covarianceModel_;
  util::Timer timer(timername_, "Constructor");
  randomizationSize_ = conf.getInt("randomization size", 0);
  fullInverse_ = conf.getBool("full inverse", false);
  fullInverseIterations_ = conf.getInt("full inverse iterations", 10);
  fullInverseAccuracy_ = conf.getDouble("full inverse accuracy", 1.0e-3);
  const eckit::LocalConfiguration cvconf = conf.getSubConfiguration("linear variable change");
  eckit::LocalConfiguration invcvconf =
        conf.getSubConfiguration("inverse linear variable change");
  if (!cvconf.empty()) {
    if (cvconf.has("input variables")) {
      BVars_.reset(new Variables(cvconf, "input variables"));
    }
    if (cvconf.has("output variables")) {
      anaVars_.reset(new Variables(cvconf, "output variables"));
    } else {
      anaVars_.reset(new Variables());
    }
    for (size_t jt = 0; jt < fg.size(); ++jt) {
      linVarChg_.push_back(std::make_unique<LinearVariableChange_>(resol, cvconf));
      linVarChg_[jt]->changeVarTraj(fg[jt], *anaVars_);
    }
    // If there is an inverse linear variable change configuration, set up
    // inverses (for use in inverseMultiply)
    if (!invcvconf.empty()) {
      // Overwrite the inverse linear variable change variables configuration
      // silently to the correct configuration.
      // TODO(someone): Remove this once model interfaces are handling
      // variable changes consistently without the need for input/output variables
      // specified in yaml.
      invcvconf.set("input variables", anaVars_->variables());
      invcvconf.set("output variables", BVars_->variables());
      for (size_t jt = 0; jt < fg.size(); ++jt) {
        invLinVarChg_.push_back(std::make_unique<LinearVariableChange_>(resol, invcvconf));
        invLinVarChg_[jt]->changeVarTraj(fg[jt], *BVars_);
      }
    }
  }
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::ModelSpaceCovarianceBase done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovarianceBase<MODEL>::randomize(Increment4D_ & dx) const {
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::randomize starting" << std::endl;
  util::Timer timer(timername_, "randomize");
  // TODO(notguillaume): Generalize to non-square change of variable
  this->doRandomize(dx);
  if (linVarChg_.size() > 0) {
    ASSERT(linVarChg_.size() == dx.size());
    for (size_t jt = 0; jt < dx.size(); ++jt) {
      linVarChg_[jt]->changeVarTL(dx[jt], *anaVars_);
    }
  }
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovarianceBase<MODEL>::multiply(const Increment4D_ & dxi,
                                               Increment4D_ & dxo) const {
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::multiply starting" << std::endl;
  util::Timer timer(timername_, "multiply");
  if (linVarChg_.size() > 0) {
    ASSERT(linVarChg_.size() == dxi.size());
    // Copy input increment and apply adjoint variable change (to control variables)
    Increment4D_ dxiTemp(dxi);
    for (size_t jt = 0; jt < dxiTemp.size(); ++jt) {
      linVarChg_[jt]->changeVarAD(dxiTemp[jt], *BVars_);
    }

    // Create temporary output increment
    Increment4D_ dxoTemp(dxiTemp, false);

    // Apply background error model
    this->doMultiply(dxiTemp, dxoTemp);

    // Apply control to analysis/model variable change
    for (size_t jt = 0; jt < dxoTemp.size(); ++jt) {
      linVarChg_[jt]->changeVarTL(dxoTemp[jt], *anaVars_);
    }
    // Copy to output increment
    dxo = dxoTemp;
  } else {
    this->doMultiply(dxi, dxo);
  }
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovarianceBase<MODEL>::inverseMultiply(const Increment4D_ & dxi,
                                                      Increment4D_ & dxo) const {
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::inverseMultiply starting" << std::endl;
  util::Timer timer(timername_, "inverseMultiply");
  if (fullInverse_) {
    // Approximate full inverse using GMRESR
    IdentityMatrix<Increment4D_> Id;
    dxo.zero();
    GMRESR(dxo, dxi, *this, Id, fullInverseIterations_, fullInverseAccuracy_);
  } else {
    if (linVarChg_.size() > 0) {
      ASSERT(linVarChg_.size() == dxi.size());
      if (invLinVarChg_.size() != linVarChg_.size()) {
        Log::error() << "ModelSpaceCovarianceBase::inverseMultiply: linear variable "
                     << "change specified, but not an inverse linear variable change."
                     << std::endl << "Specify inverse linear variable change or use "
                     << "'full inverse: true' to compute the inverse." << std::endl;
        throw eckit::BadParameter("ModelSpaceCovarianceBase::inverseMultiply failed: "
                                  "missing inverse linear variable change", Here());
      }
      // Copy input increment and apply inverse variable change (K^{-1})
      Increment4D_ dxiTemp(dxi);
      for (size_t jt = 0; jt < dxiTemp.size(); ++jt) {
        invLinVarChg_[jt]->changeVarTL(dxiTemp[jt], *anaVars_);
      }

      // Create temporary output increment
      Increment4D_ dxoTemp(dxiTemp, false);

      // Apply background error model
      this->doInverseMultiply(dxiTemp, dxoTemp);

      // Apply adjoint inverse variable change (K^T^{-1})
      for (size_t jt = 0; jt < dxoTemp.size(); ++jt) {
        invLinVarChg_[jt]->changeVarAD(dxoTemp[jt], *BVars_);
      }

      // Copy to output increment
      dxo = dxoTemp;
    } else {
      this->doInverseMultiply(dxi, dxo);
    }
  }
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovarianceBase<MODEL>::getVariance(Increment4D_ & variance) const {
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::getVariance starting" << std::endl;
  util::Timer timer(timername_, "getVariance");
  Increment_ dx(variance[0]);
  Increment_ dxsq(variance[0]);
  Increment_ mean(variance[0]);
  mean.zero();
  variance.zero();
  for (size_t ie = 0; ie < randomizationSize_; ++ie) {
    this->randomize(dx);
    dx -= mean;
    dxsq = dx;
    dxsq.schur_product_with(dx);
    double rk_var = static_cast<double>(ie)/static_cast<double>(ie+1);
    double rk_mean = 1.0/static_cast<double>(ie+1);
    variance[0].axpy(rk_var, dxsq, false);
    mean.axpy(rk_mean, dx, false);
  }
  double rk_norm = 1.0/static_cast<double>(randomizationSize_-1);
  variance *= rk_norm;
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::getVariance done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_MODELSPACECOVARIANCEBASE_H_
