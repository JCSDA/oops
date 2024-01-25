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
#include "oops/base/ModelSpaceCovarianceParametersBase.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Random.h"
#include "oops/util/Timer.h"

namespace util {
  class DateTime;
}

namespace oops {

// -----------------------------------------------------------------------------

// Should this be one with the ErrorCovariance class in the interface directory? YT

/// Abstract base class for model space error covariances.
///
/// Note: subclasses can opt to extract their settings either from a Configuration object or from a
/// subclass of Parameters.
///
/// In the former case, they should provide a constructor taking a const reference to an
/// eckit::Configuration object. In the latter case, the implementer should first define a subclass
/// of Parameters holding the settings of the model-space covariance change in question. The latter
/// should then typedef `Parameters_` to the name of that subclass and provide a constructor taking
/// a const reference to an instance of that subclass.
template <typename MODEL>
class ModelSpaceCovarianceBase {
  typedef Geometry<MODEL>                                       Geometry_;
  typedef GeometryIterator<MODEL>                               GeometryIterator_;
  typedef State4D<MODEL>                                        State4D_;
  typedef Increment<MODEL>                                      Increment_;
  typedef Increment4D<MODEL>                                    Increment4D_;
  typedef LinearVariableChange<MODEL>                           LinearVariableChange_;

 public:
  ModelSpaceCovarianceBase(const Geometry_ &, const ModelSpaceCovarianceParametersBase<MODEL> &,
                           const State4D_ &, const State4D_ &);
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
  std::unique_ptr<Variables> anaVars_;
  std::string timername_;
};

// =============================================================================

template <typename MODEL>
class CovarianceFactory;

// -----------------------------------------------------------------------------

/// \brief A subclass of ModelSpaceCovarianceParametersBase storing the values of all options in a
/// single Configuration object.
///
/// This object can be accessed by calling the value() method of the \p config member variable.
///
/// The ConfigurationParameter class does not perform any parameter validation; models using
/// GenericModelSpaceCovarianceParameters should therefore ideally be refactored, replacing this
/// class with a dedicated subclass of ModelSpaceCovarianceParametersBase storing each parameter in
/// a separate (Optional/Required)Parameter object.
///
template <typename MODEL>
class GenericModelSpaceCovarianceParameters : public ModelSpaceCovarianceParametersBase<MODEL> {
  OOPS_CONCRETE_PARAMETERS(GenericModelSpaceCovarianceParameters,
                           ModelSpaceCovarianceParametersBase<MODEL>)
 public:
  ConfigurationParameter config{this};
};

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// ModelSpaceCovarianceParametersBase.
template <typename MODEL>
class ModelSpaceCovarianceParametersWrapper : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelSpaceCovarianceParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of ModelSpaceCovarianceParametersBase
  /// controlling the behavior of a covariance model. The type of the subclass is determined by
  /// the value of the "name" key in the Configuration object from which this object is
  /// deserialized.
  RequiredPolymorphicParameter<ModelSpaceCovarianceParametersBase<MODEL>, CovarianceFactory<MODEL>>
    covarianceParameters{"covariance model", this};
};

// =============================================================================

/// Covariance Factory
template <typename MODEL>
class CovarianceFactory {
  typedef Geometry<MODEL> Geometry_;
  typedef State4D<MODEL>  State4D_;

 public:
  /// \brief Create and return a new covariance model.
  ///
  /// The covariance model is determined by the \c covarianceModel attribute of \p parameters.
  /// \p parameters must be an instance of the subclass of ModelSpaceCovarianceParametersBase
  /// associated with that covariance model, otherwise an exception will be thrown.
  static ModelSpaceCovarianceBase<MODEL> * create(const Geometry_ &,
                                                  const Variables &,
                                                  const ModelSpaceCovarianceParametersBase<MODEL> &,
                                                  const State4D_ &,
                                                  const State4D_ &);

  /// \brief Create and return a new covariance model.
  ///
  /// Deprecated overload taking a Configuration instead of a ModelSpaceCovarianceParametersBase.
  static ModelSpaceCovarianceBase<MODEL> * create(const Geometry_ &,
                                                  const Variables &,
                                                  const eckit::Configuration &,
                                                  const State4D_ &,
                                                  const State4D_ &);

  /// \brief Create and return an instance of the subclass of ModelSpaceCovarianceParametersBase
  /// storing parameters of the specified covariance model.
  static std::unique_ptr<ModelSpaceCovarianceParametersBase<MODEL>> createParameters(
      const std::string &covarianceModel);

  /// \brief Return the names of all covariance models that can be created by one of the
  /// registered makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~CovarianceFactory() = default;

 protected:
  /// \brief Register a maker able to create covariance models of type \p name.
  explicit CovarianceFactory(const std::string &name);

 private:
  virtual ModelSpaceCovarianceBase<MODEL> * make(const Geometry_ &, const Variables &,
                                                 const ModelSpaceCovarianceParametersBase<MODEL> &,
                                                 const State4D_ &, const State4D_ &) = 0;

  virtual std::unique_ptr<ModelSpaceCovarianceParametersBase<MODEL>> makeParameters() const = 0;

  static std::map < std::string, CovarianceFactory<MODEL> * > & getMakers() {
    static std::map < std::string, CovarianceFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class COVAR>
class CovarMaker : public CovarianceFactory<MODEL> {
  /// Defined as T::Parameters_ if COVAR defines a Parameters_ type; otherwise as
  /// GenericModelSpaceCovarianceParameters<MODEL>.
  typedef TParameters_IfAvailableElseFallbackType_t<
    COVAR, GenericModelSpaceCovarianceParameters<MODEL>> Parameters_;

  typedef Geometry<MODEL> Geometry_;
  typedef State4D<MODEL>  State4D_;

  ModelSpaceCovarianceBase<MODEL> * make(const Geometry_ & resol, const Variables & vars,
                                const ModelSpaceCovarianceParametersBase<MODEL> & params,
                                const State4D_ & xb, const State4D_ & fg) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new COVAR(resol, vars,
                     parametersOrConfiguration<HasParameters_<COVAR>::value>(stronglyTypedParams),
                     xb, fg);
  }

  std::unique_ptr<ModelSpaceCovarianceParametersBase<MODEL>> makeParameters() const override {
    return boost::make_unique<Parameters_>();
  }

 public:
  explicit CovarMaker(const std::string & name) : CovarianceFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
CovarianceFactory<MODEL>::CovarianceFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in covariance factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>* CovarianceFactory<MODEL>::create(
    const Geometry_ & resol,
    const Variables & vars,
    const ModelSpaceCovarianceParametersBase<MODEL> & parameters,
    const State4D_ & xb,
    const State4D_ & fg) {
  const std::string id = parameters.covarianceModel.value().value();
  Log::trace() << "ModelSpaceCovarianceBase type = " << id << std::endl;
  typename std::map<std::string, CovarianceFactory<MODEL>*>::iterator jcov = getMakers().find(id);
  if (jcov == getMakers().end()) {
    Log::error() << id << " does not exist in CovarianceFactory." << std::endl;
    Log::error() << "CovarianceFactory has " << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, CovarianceFactory<MODEL>*>::const_iterator
         jj = getMakers().begin(); jj != getMakers().end(); ++jj) {
       Log::error() << "A " << jj->first << " B" << std::endl;
    }
    throw std::runtime_error(id + " does not exist in covariance factory.");
  }
  Variables vars_in{vars};
  if (parameters.variableChange.value() != boost::none) {
    if (parameters.variableChange.value()->inputVariables.value() != boost::none) {
        vars_in = *(parameters.variableChange.value()->inputVariables.value());
    }
  }
  return (*jcov).second->make(resol, vars_in, parameters, xb, fg);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>* CovarianceFactory<MODEL>::create(const Geometry_ & resol,
                                const Variables & vars, const eckit::Configuration & conf,
                                const State4D_ & xb, const State4D_ & fg) {
  ModelSpaceCovarianceParametersWrapper<MODEL> parameters;
  parameters.validateAndDeserialize(conf);
  return create(resol, vars, parameters.covarianceParameters, xb, fg);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<ModelSpaceCovarianceParametersBase<MODEL>>
  CovarianceFactory<MODEL>::createParameters(const std::string &name) {
  typename std::map<std::string, CovarianceFactory<MODEL>*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in CovarianceFactory");
  }
  return it->second->makeParameters();
}

// =============================================================================

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>::ModelSpaceCovarianceBase(const Geometry_ & resol,
                     const ModelSpaceCovarianceParametersBase<MODEL> & parameters,
                                                          const State4D_ & xb, const State4D_ & fg)
{
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::ModelSpaceCovarianceBase starting" << std::endl;
  const boost::optional<std::string> &covarianceModel = parameters.covarianceModel.value();
  if (covarianceModel == boost::none) {
    covarianceModel_ = "none";
  } else {
    covarianceModel_ = *covarianceModel;
  }
  timername_ = "oops::Covariance::" + covarianceModel_;
  util::Timer timer(timername_, "Constructor");
  randomizationSize_ = parameters.randomizationSize.value().get_value_or(0);
  fullInverse_ = parameters.fullInverse;
  fullInverseIterations_ = parameters.fullInverseIterations;
  fullInverseAccuracy_ = parameters.fullInverseAccuracy;
  const boost::optional<typename LinearVariableChange<MODEL>::Parameters_> &variableChangeParms =
      parameters.variableChange;
  if (variableChangeParms != boost::none) {
      const boost::optional<Variables> &inputVars = variableChangeParms->inputVariables;
      if (inputVars != boost::none) {
          BVars_.reset(new Variables(*inputVars));
      }
      const boost::optional<Variables> &outputVars = variableChangeParms->outputVariables;
      if (outputVars != boost::none) {
          anaVars_.reset(new Variables(*outputVars));
      } else {
          anaVars_.reset(new Variables());
      }
      for (size_t jt = 0; jt < fg.size(); ++jt) {
        linVarChg_.push_back(std::make_unique<LinearVariableChange_>(resol, *variableChangeParms));
        linVarChg_[jt]->changeVarTraj(fg[jt], *anaVars_);
      }
  }
  Log::trace() << "ModelSpaceCovarianceBase<MODEL>::ModelSpaceCovarianceBase done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>::ModelSpaceCovarianceBase(const Geometry_ & resol,
                                                          const eckit::Configuration & conf,
                                                          const State4D_ & xb, const State4D_ & fg)
  : ModelSpaceCovarianceBase(
      resol, validateAndDeserialize<GenericModelSpaceCovarianceParameters<MODEL>>(conf), xb, fg)
{}

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
      // Copy input increment and apply inverse variable change (K^{-1})
      Increment4D_ dxiTemp(dxi);
      for (size_t jt = 0; jt < dxiTemp.size(); ++jt) {
        linVarChg_[jt]->changeVarInverseTL(dxiTemp[jt], *BVars_);
      }

      // Create temporary output increment
      Increment4D_ dxoTemp(dxiTemp, false);

      // Apply background error model
      this->doInverseMultiply(dxiTemp, dxoTemp);

      // Apply adjoint inverse variable change (K^T^{-1})
      for (size_t jt = 0; jt < dxoTemp.size(); ++jt) {
        linVarChg_[jt]->changeVarInverseAD(dxoTemp[jt], *anaVars_);
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
