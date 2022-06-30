/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSLOCALIZATIONBASE_H_
#define OOPS_BASE_OBSLOCALIZATIONBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "oops/base/ObsLocalizationParametersBase.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

namespace oops {

/// Base class for observation-space localization.
/// Defines the interfaces for observation space localization.
/// Use this class as a base class for OBS- and MODEL-specific implementations.
template<typename MODEL, typename OBS>
class ObsLocalizationBase : public util::Printable,
                                   private boost::noncopyable {
  typedef typename MODEL::GeometryIterator          GeometryIterator_;
  typedef typename OBS::ObsVector                   ObsVector_;
 public:
  ObsLocalizationBase() = default;
  virtual ~ObsLocalizationBase() = default;

  /// compute obs-space localization: update \p locfactor with observation-space
  /// localization values between observations and \p point in model-space.
  /// update means that locfactor values that are passed in are multiplied by the
  /// locfactor computed inside of this function
  /// Set \p locfactor to missing value for observations that are not local.
  /// Method used in oops. Calls `computeLocalization` abstract method, and
  /// passes OBS- and MODEL-specific classes to the OBS- and MODEL-specific
  /// implementations of ObsLocalization.
  void computeLocalization(const GeometryIterator<MODEL> & point,
                           ObsVector<OBS> & locfactor) const {
    computeLocalization(point.geometryiter(), locfactor.obsvector());
  }

  /// compute obs-space localization: update \p locfactor with observation-space
  /// localization values between observations and \p point in model-space.
  /// Set \p locfactor to missing value for observations that are not local.
  virtual void computeLocalization(const GeometryIterator_ & point,
                                   ObsVector_ & locfactor) const = 0;
};

template <typename MODEL, typename OBS> class ObsLocalizationFactory;

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// ObsLocalizationParametersBase.
template <typename MODEL, typename OBS>
class ObsLocalizationParametersWrapper : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsLocalizationParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of ObsLocalizationParametersBase
  /// controlling the behavior of an obs-space localization. The type of the subclass is determined
  /// by the value of the "localization method" key in the Configuration object from which this
  /// object is deserialized.
  RequiredPolymorphicParameter<ObsLocalizationParametersBase,
                               ObsLocalizationFactory<MODEL, OBS>>
      obslocParameters{"localization method", this};
};

/// ObsLocalization Factory
template <typename MODEL, typename OBS>
class ObsLocalizationFactory {
  typedef ObsSpace<OBS>  ObsSpace_;

 public:
  static std::unique_ptr<ObsLocalizationBase<MODEL, OBS>> create(
                        const ObsLocalizationParametersBase &, const ObsSpace_ &);

  /// \brief Create and return an instance of the subclass of ObsLocalizationParametersBase
  /// storing parameters of obs-space localizations of the specified type.
  static std::unique_ptr<ObsLocalizationParametersBase> createParameters(
      const std::string &name);

  /// \brief Return the names of all obs localizations that can be created by one of the
  /// registered makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~ObsLocalizationFactory() = default;

 protected:
  /// \brief Register a maker able to create obs localizations of type \p name.
  explicit ObsLocalizationFactory(const std::string & name);

 private:
  virtual std::unique_ptr<ObsLocalizationBase<MODEL, OBS>> make(
                          const ObsLocalizationParametersBase &, const ObsSpace_ &) = 0;

  virtual std::unique_ptr<ObsLocalizationParametersBase> makeParameters() const = 0;

  static std::map < std::string, ObsLocalizationFactory<MODEL, OBS> * > & getMakers() {
    static std::map < std::string, ObsLocalizationFactory<MODEL, OBS> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class OBS, class T>
class ObsLocalizationMaker : public ObsLocalizationFactory<MODEL, OBS> {
  typedef ObsSpace<OBS>  ObsSpace_;
  typedef typename T::Parameters_ Parameters_;

  std::unique_ptr<ObsLocalizationBase<MODEL, OBS>> make(
                                     const ObsLocalizationParametersBase & params,
                                     const ObsSpace_ & obspace) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return std::make_unique<T>(stronglyTypedParams, obspace.obsspace());
  }

  std::unique_ptr<ObsLocalizationParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit ObsLocalizationMaker(const std::string & name) :
    ObsLocalizationFactory<MODEL, OBS>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
ObsLocalizationFactory<MODEL, OBS>::ObsLocalizationFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in obs localization factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<ObsLocalizationBase<MODEL, OBS>> ObsLocalizationFactory<MODEL, OBS>::create(
                const ObsLocalizationParametersBase & params, const ObsSpace_ & obspace) {
  Log::trace() << "ObsLocalizationBase<MODEL, OBS>::create starting" << std::endl;
  const std::string id = params.method;
  typename std::map<std::string, ObsLocalizationFactory<MODEL, OBS>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in obs localization factory.");
  }
  std::unique_ptr<ObsLocalizationBase<MODEL, OBS>> ptr(jloc->second->make(params, obspace));
  Log::trace() << "ObsLocalizationBase<MODEL, OBS>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::unique_ptr<ObsLocalizationParametersBase>
ObsLocalizationFactory<MODEL, OBS>::createParameters(const std::string &name) {
  const auto & it = getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in ObsLocalizationFactory");
  }
  return it->second->makeParameters();
}

}  // namespace oops

#endif  // OOPS_BASE_OBSLOCALIZATIONBASE_H_
