/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_ANALYTICINITBASE_H_
#define OOPS_GENERIC_ANALYTICINITBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "oops/interface/GeoVaLs.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"

namespace oops {

/// \brief Initializes GeoVaLs with analytic formula
template <typename OBS>
class AnalyticInitBase : private util::ObjectCounter<AnalyticInitBase<OBS> > {
  typedef GeoVaLs<OBS>                GeoVaLs_;
  typedef SampledLocations<OBS>       SampledLocations_;

 public:
  static const std::string classname() {return "oops::AnalyticInitBase";}

  AnalyticInitBase() = default;
  virtual ~AnalyticInitBase() = default;

/*! \brief Fill GeoVaLs with values computed by analytic function.
 *
 * \details **AnalyticInit::fillGeoVaLs** was introduced in May, 2018 (initially
 * as a GeoVaLs constructor) for use with the interpolation test. The interpolation test
 * requires an initialization of a GeoVaLs object based on the same analytic
 * formulae used for the State initialization (see test::TestStateInterpolation()
 * for further information).  This in turn requires information about the
 * vertical profile in addition to the latitude and longitude positional
 * information in the SampledLocations object.  Currently, this information
 * about the vertical profile is obtained from an existing GeoVaLs object
 * (passed as *gvals*) that represents the output of the State::interpolate()
 * method.
 *
 * \date May, 2018: created as a constructor (M. Miesch, JCSDA)
 * \date June, 2018: moved to a method (M. Miesch, JCSDA)
 *
 * \sa test::TestStateInterpolation()
 */
  virtual void fillGeoVaLs(const SampledLocations_ &, GeoVaLs_ &) const = 0;
};

// -----------------------------------------------------------------------------

template <typename OBS>
class AnalyticInitFactory;

// -----------------------------------------------------------------------------

/// \brief Configuration parameters of an implementation of analytic init
class AnalyticInitParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(AnalyticInitParametersBase, Parameters)
 public:
  /// \brief Name of the analytic init method
  RequiredParameter<std::string> method{"method", this};
};

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// AnalyticInitParametersBase.
template <typename OBS>
class AnalyticInitParametersWrapper : public Parameters {
  OOPS_CONCRETE_PARAMETERS(AnalyticInitParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of AnalyticInitParametersBase
  /// controlling the behavior of an analytic init. The type of the subclass is
  /// determined by the value of the "method" key in the Configuration object from which
  /// this object is deserialized.
  RequiredPolymorphicParameter<AnalyticInitParametersBase, AnalyticInitFactory<OBS>>
    analyticInitParameters{"method", this};
};

// -----------------------------------------------------------------------------

/// A factory creating instances of concrete subclasses of AnalyticInitBase.
template <typename OBS>
class AnalyticInitFactory {
  typedef AnalyticInitBase<OBS> AnalyticInitBase_;

 public:
  /// \brief Create and return a new analytic init
  ///
  /// The analytic init's type is determined by the `method` attribute of \p params.
  /// \p params must be an instance of the subclass of AnalyticInitParametersBase
  /// associated with that method, otherwise an exception will be thrown.
  static std::unique_ptr<AnalyticInitBase_> create(const AnalyticInitParametersBase &params);

  /// \brief Create and return an instance of the subclass of AnalyticInitParametersBase
  /// storing parameters of analytic init method of the specified type.
  static std::unique_ptr<AnalyticInitParametersBase> createParameters(const std::string &method);

  /// \brief Return the names of all methods that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~AnalyticInitFactory() = default;

 protected:
  /// \brief Register a maker able to create analytic init of type \p method.
  explicit AnalyticInitFactory(const std::string &method);

 private:
  virtual std::unique_ptr<AnalyticInitBase_> make(const AnalyticInitParametersBase &) = 0;

  virtual std::unique_ptr<AnalyticInitParametersBase> makeParameters() const = 0;

  static std::map < std::string, AnalyticInitFactory<OBS> * > & getMakers() {
    static std::map < std::string, AnalyticInitFactory<OBS> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

/// \brief A subclass of AnalyticInitFactory able to create instances of T (a concrete subclass of
/// AnalyticInitBase<OBS>).
template<class OBS, class T>
class AnalyticInitMaker : public AnalyticInitFactory<OBS> {
  typedef AnalyticInitBase<OBS>       AnalyticInitBase_;
  typedef typename T::Parameters_ Parameters_;

  std::unique_ptr<AnalyticInitBase_> make(const AnalyticInitParametersBase & parameters) override {
    const auto &stronglyTypedParameters = dynamic_cast<const Parameters_&>(parameters);
    return std::make_unique<T>(stronglyTypedParameters);
  }

  std::unique_ptr<AnalyticInitParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit AnalyticInitMaker(const std::string & name) : AnalyticInitFactory<OBS>(name) {}
};

// -----------------------------------------------------------------------------


template <typename OBS>
AnalyticInitFactory<OBS>::AnalyticInitFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in analytic init factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::unique_ptr<AnalyticInitBase<OBS>>
AnalyticInitFactory<OBS>::create(const AnalyticInitParametersBase & params) {
  Log::trace() << "AnalyticInitFactory<OBS>::create starting" << std::endl;
  const std::string &id = params.method;
  Log::trace() << "AnalyticInit type is: " << id << std::endl;
  typename std::map<std::string, AnalyticInitFactory<OBS>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in analytic init factory.");
  }
  std::unique_ptr<AnalyticInitBase<OBS>> ptr(jerr->second->make(params));
  Log::trace() << "AnalyticInitFactory<OBS>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::unique_ptr<AnalyticInitParametersBase> AnalyticInitFactory<OBS>::createParameters(
    const std::string &name) {
  Log::trace() << "AnalyticInitFactory<OBS>::createParameters starting" << std::endl;
  typename std::map<std::string, AnalyticInitFactory<OBS>*>::iterator it = getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in analytic init factory");
  }
  return it->second->makeParameters();
}

}  // namespace oops

#endif  // OOPS_GENERIC_ANALYTICINITBASE_H_
