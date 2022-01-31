/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_OBSFILTERBASE_H_
#define OOPS_GENERIC_OBSFILTERBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/generic/ObsFilterParametersBase.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

namespace oops {

template <typename OBS> class GeoVaLs;
template <typename OBS> class ObsDiagnostics;
template <typename OBS> class ObsSpace;
template <typename OBS> class ObsVector;
template <typename OBS, typename DATA> class ObsDataVector;

/// \brief Stage at which to run a filter.
///
/// AUTO indicates that the stage will be determined automatially based on the data requested
/// by the filter.
/// PRE, PRIOR and POST indicate that the filter will be run at the relevant stage
/// but a check will be performed to ensure that all of the data structures required by the filter
/// are present.
enum class FilterStage {AUTO,
                        PRE,
                        PRIOR,
                        POST};

/// \brief Base class for generic implementations of filters processing observations.
///
/// Use this class as a base class for generic implementations
/// and interface::ObsFilterBase as a base class for OBS-specific implementations.
///
/// Note: implementations of this interface can opt to extract their settings either from
/// a Configuration object or from a subclass of ObsFilterParametersBase.
///
/// In the former case, they should provide a constructor with the following signature:
///
///    ObsFilter(const ObsSpace_ &, const eckit::Configuration &,
///              ObsDataPtr_<int>, ObsDataPtr_<float>);
///
/// In the latter case, the implementer should first define a subclass of ObsFilterParametersBase
/// holding the settings of the filter in question. The implementation of the ObsFilter interface
/// should then typedef `Parameters_` to the name of that subclass and provide a constructor with
/// the following signature:
///
///    ObsFilter(const ObsSpace_ &, const Parameters_ &,
///              ObsDataPtr_<int>, ObsDataPtr_<float>);
template <typename OBS>
class ObsFilterBase : public util::Printable,
                      private boost::noncopyable {
  typedef GeoVaLs<OBS>             GeoVaLs_;
  typedef ObsDiagnostics<OBS>      ObsDiags_;
  typedef ObsVector<OBS>           ObsVector_;

 public:
  ObsFilterBase() {}
  virtual ~ObsFilterBase() {}

  /// \brief Perform any observation processing steps that do not require access to GeoVaLs or
  /// outputs produced by the observation operator.
  virtual void preProcess() = 0;

  /// \brief Perform any observation processing steps that require access to GeoVaLs, but not to
  /// outputs produced by the observation operator.
  virtual void priorFilter(const GeoVaLs_ &gv) = 0;

  /// \brief Perform any observation processing steps that require access to both GeoVaLs and
  /// outputs produced by the observation operator.
  ///
  /// \param gv
  ///   GeoVaLs.
  /// \param ov
  ///   Model equivalents produced by the observation operator.
  /// \param bv
  ///   Bias of departure produced by the observation operator.
  /// \param dv
  ///   Observation diagnostics produced by the observation operator.
  virtual void postFilter(const GeoVaLs_ & gv,
                          const ObsVector_ &ov,
                          const ObsVector_ &bv,
                          const ObsDiags_ &dv) = 0;

  /// \brief Check the required filter data are present prior to running this filter.
  virtual void checkFilterData(const FilterStage filterStage) = 0;

  /// \brief Return the list of GeoVaLs required by this filter.
  virtual Variables requiredVars() const = 0;

  /// \brief Return the list of observation diagnostics required by this filter.
  virtual Variables requiredHdiagnostics() const = 0;
};

// =============================================================================

template <typename OBS>
class FilterFactory;

// -----------------------------------------------------------------------------

/// \brief A subclass of ObsFilterParametersBase storing the values of all options in a
/// single Configuration object.
///
/// This object can be accessed by calling the value() method of the \p config member variable.
///
/// The ConfigurationParameter class does not perform any parameter validation; filters using
/// GenericFilterParameters should therefore ideally be refactored, replacing this
/// class with a dedicated subclass of ObsFilterParametersBase storing each parameter in
/// a separate (Optional/Required)Parameter object.
class GenericObsFilterParameters : public ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericObsFilterParameters, ObsFilterParametersBase)
 public:
  ConfigurationParameter config{this};
};

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// ObsFilterParametersBase.
template <typename OBS>
class ObsFilterParametersWrapper : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsFilterParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of ObsFilterParametersBase
  /// controlling the behavior of an observation filter. The type of the subclass is determined
  /// by the value of the "filter" key in the Configuration object from which this object
  /// is deserialized.
  RequiredPolymorphicParameter<ObsFilterParametersBase, FilterFactory<OBS>> filterParameters{
      "filter", this};

  /// Indices of iterations at which this filter should be applied.
  OptionalParameter<std::string> applyAtIterations{"apply at iterations", this};
};

// =============================================================================

/// ObsFilter factory
template <typename OBS>
class FilterFactory {
  typedef ObsSpace<OBS>    ObsSpace_;
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ObsDataVector<OBS, DATA> >;

 public:
  /// \brief Create and return a new observation filter.
  ///
  /// The type of the filter is determined by the `filter` attribute of \p params. \p params
  /// must be an instance of the subclass of ObsFilterParametersBase associated with that filter,
  /// otherwise an exception will be thrown.
  static std::unique_ptr<ObsFilterBase<OBS>> create(const ObsSpace_ &,
                                                    const ObsFilterParametersBase & params,
                                                    ObsDataPtr_<int> flags = ObsDataPtr_<int>(),
                                               ObsDataPtr_<float> obserr = ObsDataPtr_<float>());

  /// \brief Create and return an instance of the subclass of ObsFilterParametersBase
  /// storing parameters of observation filters of the specified type.
  static std::unique_ptr<ObsFilterParametersBase> createParameters(
      const std::string &name);

  /// \brief Return the names of all filters that can be created by one of the
  /// registered makers.
  static std::vector<std::string> getMakerNames() {
    return keys(getMakers());
  }

  virtual ~FilterFactory() = default;

 protected:
  /// \brief Register a maker able to create observation filters of type \p name.
  explicit FilterFactory(const std::string &name);

 private:
  virtual std::unique_ptr<ObsFilterBase<OBS>> make(const ObsSpace_ &,
                                                   const ObsFilterParametersBase &,
                                                   ObsDataPtr_<int> &, ObsDataPtr_<float> &) = 0;

  virtual std::unique_ptr<ObsFilterParametersBase> makeParameters() const = 0;

  static std::map < std::string, FilterFactory<OBS> * > & getMakers() {
    static std::map < std::string, FilterFactory<OBS> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class OBS, class T>
class FilterMaker : public FilterFactory<OBS> {
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// GenericObsFilterParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<T, GenericObsFilterParameters> Parameters_;

  typedef ObsSpace<OBS>    ObsSpace_;
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ObsDataVector<OBS, DATA> >;

  std::unique_ptr<ObsFilterBase<OBS>> make(const ObsSpace_ & os,
                                           const ObsFilterParametersBase & params,
                                           ObsDataPtr_<int> & flags,
                                           ObsDataPtr_<float> & obserr) override {
        const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
        return std::make_unique<T>(os,
                     parametersOrConfiguration<HasParameters_<T>::value>(stronglyTypedParams),
                     flags,
                     obserr);
  }

  std::unique_ptr<ObsFilterParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit FilterMaker(const std::string & name) : FilterFactory<OBS>(name) {}
};

// =============================================================================

template <typename OBS>
FilterFactory<OBS>::FilterFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in obs filter factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::unique_ptr<ObsFilterBase<OBS>>
FilterFactory<OBS>::create(const ObsSpace_ & os, const ObsFilterParametersBase & params,
                           ObsDataPtr_<int> flags, ObsDataPtr_<float> obserr) {
  Log::trace() << "FilterFactory<OBS>::create starting" << std::endl;
  const std::string &id = params.filter.value().value();
  typename std::map<std::string, FilterFactory<OBS>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in obs filter factory." << std::endl;
    Log::error() << "Obs Filter Factory has " << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, FilterFactory<OBS>*>::const_iterator
         jj = getMakers().begin(); jj != getMakers().end(); ++jj) {
       Log::error() << "A " << jj->first << " Filter" << std::endl;
    }
    throw std::runtime_error(id + " does not exist in obs filter factory.");
  }
  std::unique_ptr<ObsFilterBase<OBS>> ptr(jloc->second->make(os, params, flags, obserr));
  Log::trace() << "FilterFactory<OBS>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::unique_ptr<ObsFilterParametersBase>
FilterFactory<OBS>::createParameters(const std::string &name) {
  typename std::map<std::string, FilterFactory<OBS>*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in ObsFilterFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_OBSFILTERBASE_H_
