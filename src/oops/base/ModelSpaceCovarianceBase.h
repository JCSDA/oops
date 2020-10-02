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
#include <memory>
#include <string>
#include <vector>

#include <boost/make_unique.hpp>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/ModelSpaceCovarianceParametersBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"

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
  typedef State<MODEL>                                          State_;
  typedef Increment<MODEL>                                      Increment_;
  typedef LinearVariableChangeBase<MODEL>                       LinearVariableChangeBase_;
  typedef typename boost::ptr_vector<LinearVariableChangeBase_> ChvarVec_;
  typedef typename ChvarVec_::iterator iter_;
  typedef typename ChvarVec_::const_iterator icst_;
  typedef typename ChvarVec_::const_reverse_iterator ircst_;

 public:
  ModelSpaceCovarianceBase(const State_ &, const State_ &,
                           const Geometry_ &, const ModelSpaceCovarianceParametersBase<MODEL> &);
  ModelSpaceCovarianceBase(const State_ &, const State_ &,
                           const Geometry_ &, const eckit::Configuration &);
  virtual ~ModelSpaceCovarianceBase() {}

  void randomize(Increment_ &) const;
  void multiply(const Increment_ &, Increment_ &) const;
  void inverseMultiply(const Increment_ &, Increment_ &) const;

 private:
  virtual void doRandomize(Increment_ &) const = 0;
  virtual void doMultiply(const Increment_ &, Increment_ &) const = 0;
  virtual void doInverseMultiply(const Increment_ &, Increment_ &) const = 0;

  ChvarVec_ chvars_;
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
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;

 public:
  /// \brief Create and return a new covariance model.
  ///
  /// The covariance model is determined by the \c covarianceModel attribute of \p parameters.
  /// \p parameters must be an instance of the subclass of ModelSpaceCovarianceParametersBase
  /// associated with that covariance model, otherwise an exception will be thrown.
  static ModelSpaceCovarianceBase<MODEL> * create(const ModelSpaceCovarianceParametersBase<MODEL> &,
                                                  const Geometry_ &, const Variables &,
                                                  const State_ &, const State_ &);

  /// \brief Create and return a new covariance model.
  ///
  /// Deprecated overload taking a Configuration instead of a ModelSpaceCovarianceParametersBase.
  static ModelSpaceCovarianceBase<MODEL> * create(const eckit::Configuration &,
                                                  const Geometry_ &, const Variables &,
                                                  const State_ &, const State_ &);

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
  virtual ModelSpaceCovarianceBase<MODEL> * make(const ModelSpaceCovarianceParametersBase<MODEL> &,
                                                 const Geometry_ &, const Variables &,
                                                 const State_ &, const State_ &) = 0;

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

  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;

  ModelSpaceCovarianceBase<MODEL> * make(
      const ModelSpaceCovarianceParametersBase<MODEL> & params,
      const Geometry_ & resol, const Variables & vars,
      const State_ & xb, const State_ & fg) override {
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
    Log::error() << name << " already registered in covariance factory." << std::endl;
    ABORT("Element already registered in CovarianceFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>* CovarianceFactory<MODEL>::create(
    const ModelSpaceCovarianceParametersBase<MODEL> & parameters,
    const Geometry_ & resol,
    const Variables & vars,
    const State_ & xb, const State_ & fg) {
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
    ABORT("Element does not exist in CovarianceFactory.");
  }
  Variables vars_in(vars);
  Variables vars_out;
  for (const LinearVariableChangeParametersWrapper<MODEL> &variableChange :
       boost::adaptors::reverse(parameters.variableChanges.value())) {
    const LinearVariableChangeParametersBase &variableChangeParameters =
      variableChange.variableChangeParameters;
    if (variableChangeParameters.inputVariables.value() != boost::none &&
        variableChangeParameters.outputVariables.value() != boost::none) {
      vars_out = *variableChangeParameters.outputVariables.value();
      if (!(vars_in == vars_out)) {
        Log::error() << "Input variables:  " << vars_in << std::endl;
        Log::error() << "Output variables: " << vars_out << std::endl;
        ABORT("Sequence of variable changes is not consistent");
      }
      vars_in = *variableChangeParameters.inputVariables.value();
    }
  }
  return (*jcov).second->make(parameters, resol, vars_in, xb, fg);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>* CovarianceFactory<MODEL>::create(
    const eckit::Configuration & conf,
    const Geometry_ & resol,
    const Variables & vars,
    const State_ & xb, const State_ & fg) {
  ModelSpaceCovarianceParametersWrapper<MODEL> parameters;
  parameters.validateAndDeserialize(conf);
  return create(parameters.covarianceParameters, resol, vars, xb, fg);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<ModelSpaceCovarianceParametersBase<MODEL>>
  CovarianceFactory<MODEL>::createParameters(
    const std::string &name) {
  typename std::map<std::string, CovarianceFactory<MODEL>*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in CovarianceFactory");
  }
  return it->second->makeParameters();
}

// =============================================================================

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>::ModelSpaceCovarianceBase(
    const State_ & bg, const State_ & fg,
    const Geometry_ & resol,
    const ModelSpaceCovarianceParametersBase<MODEL> & parameters) {
  for (const LinearVariableChangeParametersWrapper<MODEL> &variableChange :
       parameters.variableChanges.value()) {
    chvars_.push_back(LinearVariableChangeFactory<MODEL>::create(
                        bg, fg, resol, variableChange.variableChangeParameters));
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelSpaceCovarianceBase<MODEL>::ModelSpaceCovarianceBase(
    const State_ & bg, const State_ & fg,
    const Geometry_ & resol,
    const eckit::Configuration & conf)
  : ModelSpaceCovarianceBase(
      bg, fg, resol, validateAndDeserialize<GenericModelSpaceCovarianceParameters<MODEL>>(conf))
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovarianceBase<MODEL>::randomize(Increment_ & dx) const {
  // TODO(notguillaume): Generalize to non-square change of variable
  if (chvars_.size()) {
    this->doRandomize(dx);   // dx = C^1/2 dx
    // K_N K_N-1 ... K_1
    for (icst_ it = chvars_.begin(); it != chvars_.end(); ++it) {
      dx = it->multiply(dx);  // dx = K_i dx
    }
  } else {
    this->doRandomize(dx);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelSpaceCovarianceBase<MODEL>::multiply(const Increment_ & dxi,
                                               Increment_ & dxo) const {
  if (chvars_.size()) {
    // K_1^T K_2^T .. K_N^T
    std::unique_ptr<Increment_> dxchvarin(new Increment_(dxi));
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
    std::unique_ptr<Increment_> dxchvarin(new Increment_(dxi));
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
