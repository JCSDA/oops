/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_GENERIC_OBSERRORBASE_H_
#define OOPS_GENERIC_OBSERRORBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/ObsVector.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/PolymorphicParameter.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Base class for generic implementations of observation error covariance matrices.
///
/// Use this class as a base class for generic implementations,
/// and interface::ObsErrorBase as a base class for OBS-specific implementations.
///
/// Note: each generic implementation should typedef `Parameters_` to the name of a subclass of
/// ObsErrorParametersBase holding its configuration settings and provide a constructor with the
/// following signature:
///
///     ObsErrorBase(const Parameters_ &, ObsSpace<OBS> &);
template<typename OBS>
class ObsErrorBase : public util::Printable,
                     private boost::noncopyable {
  typedef ObsVector<OBS>        ObsVector_;
  typedef ObsSpace<OBS>         ObsSpace_;

 public:
  ObsErrorBase() = default;
  virtual ~ObsErrorBase() = default;

/// Multiply a Departure \p dy by \f$R\f$.
  virtual void multiply(ObsVector_ & dy) const = 0;
/// Multiply a Departure \p dy by \f$R^{-1}\f$.
  virtual void inverseMultiply(ObsVector_ & dy) const = 0;

/// Generate random perturbation in \p dy.
  virtual void randomize(ObsVector_ & dy) const = 0;

/// Save obs errors to the group \p name.
  virtual void save(const std::string &name) const = 0;

/// Return a copy of obs error std. dev. If this ObsVector_ is modified (e.g. by obs filters),
/// it should be passed back to update() to ensure the covariance matrix stays consistent.
  virtual ObsVector_ obserrors() const = 0;

/// Set the diagonal of the covariance matrix to \p stddev squared.
  virtual void update(const ObsVector_ &stddev) = 0;

/// Return the vector of inverse obs error variances.
  virtual ObsVector_ inverseVariance() const = 0;

/// Get mean error for Jo table.
  virtual double getRMSE() const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// =============================================================================

template <typename OBS>
class ObsErrorFactory;

// -----------------------------------------------------------------------------

/// \brief Configuration parameters of an implementation of an observation error covariance matrix
/// model.
class ObsErrorParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(ObsErrorParametersBase, Parameters)
 public:
  /// \brief Name of the covariance model.
  Parameter<std::string> model{"covariance model", "diagonal", this};
};

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// ObsErrorParametersBase.
template <typename OBS>
class ObsErrorParametersWrapper : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of ObsErrorParametersBase controlling
  /// the behavior of an observation error covariance matrix model. The type of the subclass is
  /// determined by the value of the "covariance model" key in the Configuration object from which
  /// this object is deserialized.
  PolymorphicParameter<ObsErrorParametersBase, ObsErrorFactory<OBS>>
    obsErrorParameters{"covariance model", "diagonal", this};
};

// -----------------------------------------------------------------------------

/// A factory creating instances of concrete subclasses of ObsErrorBase.
template <typename OBS>
class ObsErrorFactory {
  typedef ObsErrorBase<OBS> ObsErrorBase_;
  typedef ObsSpace<OBS>     ObsSpace_;

 public:
  /// \brief Create and return a new observation error covariance matrix model.
  ///
  /// The model's type is determined by the `covariance model` attribute of \p params.
  /// \p params must be an instance of the subclass of ObsErrorParametersBase
  /// associated with that model type, otherwise an exception will be thrown.
  static std::unique_ptr<ObsErrorBase_> create(const ObsErrorParametersBase &params,
                                               const ObsSpace_ &);

  /// \brief Create and return an instance of the subclass of ObsErrorParametersBase
  /// storing parameters of covariance matrix models of the specified type.
  static std::unique_ptr<ObsErrorParametersBase> createParameters(const std::string &model);

  /// \brief Return the names of all models that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~ObsErrorFactory() = default;

 protected:
  /// \brief Register a maker able to create covariance matrix models of type \p model.
  explicit ObsErrorFactory(const std::string &model);

 private:
  virtual std::unique_ptr<ObsErrorBase_> make(const ObsErrorParametersBase &,
                                              const ObsSpace_ &) = 0;

  virtual std::unique_ptr<ObsErrorParametersBase> makeParameters() const = 0;

  static std::map < std::string, ObsErrorFactory<OBS> * > & getMakers() {
    static std::map < std::string, ObsErrorFactory<OBS> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

/// \brief A subclass of ObsErrorFactory able to create instances of T (a concrete subclass of
/// ObsErrorBase<OBS>). Passes ObsSpace<OBS> to the constructor of T.
template<class OBS, class T>
class ObsErrorMaker : public ObsErrorFactory<OBS> {
  typedef ObsErrorBase<OBS>       ObsErrorBase_;
  typedef ObsSpace<OBS>           ObsSpace_;
  typedef typename T::Parameters_ Parameters_;

  std::unique_ptr<ObsErrorBase_> make(const ObsErrorParametersBase & parameters,
                                      const ObsSpace_ & obs) override {
    const auto &stronglyTypedParameters = dynamic_cast<const Parameters_&>(parameters);
    return std::make_unique<T>(stronglyTypedParameters, obs);
  }

  std::unique_ptr<ObsErrorParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit ObsErrorMaker(const std::string & name) : ObsErrorFactory<OBS>(name) {}
};

// =============================================================================

template <typename OBS>
ObsErrorFactory<OBS>::ObsErrorFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in obs error factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::unique_ptr<ObsErrorBase<OBS>>
ObsErrorFactory<OBS>::create(const ObsErrorParametersBase & params, const ObsSpace_ & obs) {
  Log::trace() << "ObsErrorFactory<OBS>::create starting" << std::endl;
  const std::string &id = params.model;
  Log::trace() << "ObsError matrix type is: " << id << std::endl;
  typename std::map<std::string, ObsErrorFactory<OBS>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in obs error factory.");
  }
  std::unique_ptr<ObsErrorBase<OBS>> ptr(jerr->second->make(params, obs));
  Log::trace() << "ObsErrorFactory<OBS>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::unique_ptr<ObsErrorParametersBase> ObsErrorFactory<OBS>::createParameters(
    const std::string &name) {
  Log::trace() << "ObsErrorFactory<OBS>::createParameters starting" << std::endl;
  typename std::map<std::string, ObsErrorFactory<OBS>*>::iterator it = getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in obs error factory");
  }
  std::unique_ptr<ObsErrorParametersBase> ptr(it->second->makeParameters());
  Log::trace() << "ObsErrorFactory<OBS>::createParameters done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_OBSERRORBASE_H_
