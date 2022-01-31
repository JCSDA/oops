/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSFILTERBASE_H_
#define OOPS_INTERFACE_OBSFILTERBASE_H_

#include <memory>
#include <string>

#include "oops/base/ObsVector.h"
#include "oops/generic/ObsFilterBase.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsSpace.h"

namespace oops {

namespace interface {

// -----------------------------------------------------------------------------

/// \brief Base class for OBS-specific implementations of the ObsFilter interface.
/// interface::ObsFilterBase overrides oops::ObsFilterBase methods to pass OBS-specific
/// implementations of GeoVaLs, ObsDiags, ObsSpace and ObsVector to the OBS-specific
/// implementation of ObsFilterBase.
///
/// Note: implementations of the ObsFilter interface can opt to extract their settings either from
/// a Configuration object or from a subclass of ObsFilterParametersBase.
///
/// In the former case, they should provide a constructor with the following signature:
///
///    ObsFilter(const ObsSpace_ &,  const eckit::Configuration &,
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
class ObsFilterBase : public oops::ObsFilterBase<OBS> {
  typedef typename OBS::GeoVaLs                        GeoVaLs_;
  typedef typename OBS::ObsDiagnostics                 ObsDiags_;
  typedef typename OBS::ObsSpace                       ObsSpace_;
  typedef typename OBS::ObsVector                      ObsVector_;

 public:
  // Overrides of oops::ObsFilterBase methods, converting between oops interface classes and
  // OBS-specific classes operated upon by OBS-specific implementations of the ObsFilter interface.

  void priorFilter(const GeoVaLs<OBS> &gv) final {
    this->priorFilter(gv.geovals());
  }

  void postFilter(const GeoVaLs<OBS> & gv,
                  const oops::ObsVector<OBS> &ov,
                  const oops::ObsVector<OBS> &bv,
                  const ObsDiagnostics<OBS> &dv) final {
    this->postFilter(gv.geovals(),
                     ov.obsvector(),
                     bv.obsvector(),
                     dv.obsdiagnostics());
  }

  // The methods below need to be overridden in subclasses (along with preProcess(), requiredVars()
  // and requiredHdiagnostics() and print(), methods inherited from parent classes that neither
  // take nor return instances of oops interface classes and therefore remain abstract).

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
  /// \param bias
  ///   Bias of departure produced by the observation operator.
  /// \param dv
  ///   Observation diagnostics produced by the observation operator.
  virtual void postFilter(const GeoVaLs_ & gv,
                          const ObsVector_ &ov,
                          const ObsVector_ &bv,
                          const ObsDiags_ &dv) = 0;

  /// \brief Check the required filter data are present prior to running this filter.
  virtual void checkFilterData(const oops::FilterStage filterStage) = 0;
};

// -----------------------------------------------------------------------------

/// \brief A subclass of FilterFactory able to create instances of T (a concrete subclass of
/// interface::ObsFilterBase<OBS>). Passes OBS::ObsSpace to the constructor of T.
template <typename OBS, typename T>
class FilterMaker : public FilterFactory<OBS> {
 private:
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// GenericObsFilterParameters.
  typedef TParameters_IfAvailableElseFallbackType_t<T, GenericObsFilterParameters> Parameters_;

 public:
  typedef oops::ObsSpace<OBS>   ObsSpace_;
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ObsDataVector<OBS, DATA> >;

  explicit FilterMaker(const std::string & name) : FilterFactory<OBS>(name) {}

  std::unique_ptr<oops::ObsFilterBase<OBS>> make(const ObsSpace_ & os,
                                                 const ObsFilterParametersBase & params,
                                                 ObsDataPtr_<int> & flags,
                                                 ObsDataPtr_<float> & obserr) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    const auto &paramsOrConfig =
        parametersOrConfiguration<HasParameters_<T>::value>(stronglyTypedParams);
    return std::make_unique<T>(os.obsspace(),
                               paramsOrConfig,
                               flags ? flags->obsdatavectorptr() : nullptr,
                               obserr ? obserr->obsdatavectorptr() : nullptr);
  }

  std::unique_ptr<ObsFilterParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }
};

}  // namespace interface

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSFILTERBASE_H_
