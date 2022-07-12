/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSFILTER_H_
#define OOPS_BASE_OBSFILTER_H_

#include <memory>
#include <string>

#include "oops/base/Variables.h"
#include "oops/generic/ObsFilterBase.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

template <typename OBS> class GeoVaLs;
template <typename OBS> class ObsDiagnostics;
template <typename OBS> class ObsSpace;
template <typename OBS> class ObsVector;
template <typename OBS, typename DATA> class ObsDataVector;

// -----------------------------------------------------------------------------

/// \brief A filter processing observations.
///
/// Note: to see methods that need to be implemented in a generic observation filter
/// implementation, see the ObsFilterBase class in generic/ObsFilterBase.h. To see methods that
/// need to be implemented in an OBS-specific observation filter implementation, see the
/// interface::ObsFilterBase class in interface/ObsFilterBase.h.
template <typename OBS>
class ObsFilter : public util::Printable,
                  private util::ObjectCounter<ObsFilter<OBS> > {
  typedef GeoVaLs<OBS>             GeoVaLs_;
  typedef ObsDiagnostics<OBS>      ObsDiags_;
  typedef ObsFilterBase<OBS>       ObsFilterBase_;
  typedef ObsSpace<OBS>            ObsSpace_;
  typedef ObsVector<OBS>           ObsVector_;
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ObsDataVector<OBS, DATA> >;

 public:
  static const std::string classname() {return "oops::ObsFilter";}

  /// \brief Create a new observation filter.
  ///
  /// \param obsspace
  ///   Space containing the observations to process.
  /// \param params
  ///   The filter's configuration parameters.
  /// \param qcflags
  ///   Quality control flags. They may be modified by the filter.
  /// \param obserrors
  ///   Estimates of the standard deviations of observation errors. They may be modified by the
  ///   filter.
  ObsFilter(const ObsSpace_ &obsspace, const ObsFilterParametersBase &params,
            ObsDataPtr_<int> qcflags, ObsDataPtr_<float> obserrors);
  ObsFilter(const ObsFilter &) = delete;
  ObsFilter(ObsFilter &&) = default;
  ObsFilter& operator=(const ObsFilter &) = delete;
  ObsFilter& operator=(ObsFilter &&) = default;
  ~ObsFilter();

  /// \brief Perform any observation processing steps that do not require access to GeoVaLs or
  /// outputs produced by the observation operator.
  void preProcess();

  /// \brief Perform any observation processing steps that require access to GeoVaLs, but not to
  /// outputs produced by the observation operator.
  void priorFilter(const GeoVaLs_ &);

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
  void postFilter(const GeoVaLs_ & gv,
                  const ObsVector_ &ov,
                  const ObsVector_ &bv,
                  const ObsDiags_ &dv);

  /// \brief Check the required filter data are present prior to running this filter.
  void checkFilterData(const FilterStage filterStage);

  /// \brief Return the list of GeoVaLs required by this filter.
  Variables requiredVars() const;

  /// \brief Return the list of observation diagnostics required by this filter.
  Variables requiredHdiagnostics() const;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<ObsFilterBase_> ofilt_;
  std::string filterName_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsFilter<OBS>::ObsFilter(const ObsSpace_ & os,
                          const ObsFilterParametersBase & parameters,
                          ObsDataPtr_<int> flags, ObsDataPtr_<float> obserr)
  : ofilt_(), filterName_("oops::ObsFilter::"+parameters.filter.value().value())
{
  Log::trace() << "ObsFilter<OBS>::ObsFilter starting" << std::endl;
  util::Timer timer(classname(), "ObsFilter");
  util::Timer timef(filterName_, "ObsFilter");
  ofilt_ = FilterFactory<OBS>::create(os, parameters, flags, obserr);
  Log::trace() << "ObsFilter<OBS>::ObsFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsFilter<OBS>::~ObsFilter() {
  Log::trace() << "ObsFilter<OBS>::~ObsFilter starting" << std::endl;
  util::Timer timer(classname(), "~ObsFilter");
  ofilt_.reset();
  Log::trace() << "ObsFilter<OBS>::~ObsFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsFilter<OBS>::preProcess() {
  Log::trace() << "ObsFilter<OBS>::preProcess starting" << std::endl;
  util::Timer timer(classname(), "preProcess");
  util::Timer timef(filterName_, "preProcess");
  ofilt_->preProcess();
  Log::trace() << "ObsFilter<OBS>::preProcess done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsFilter<OBS>::priorFilter(const GeoVaLs_ & gv) {
  Log::trace() << "ObsFilter<OBS>::priorFilter starting" << std::endl;
  util::Timer timer(classname(), "priorFilter");
  util::Timer timef(filterName_, "priorFilter");
  ofilt_->priorFilter(gv);
  Log::trace() << "ObsFilter<OBS>::priorFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsFilter<OBS>::postFilter(const GeoVaLs_ & gv,
                                const ObsVector_ & ov,
                                const ObsVector_ & bv,
                                const ObsDiags_ & dv) {
  Log::trace() << "ObsFilter<OBS>::postFilter starting" << std::endl;
  util::Timer timer(classname(), "postFilter");
  util::Timer timef(filterName_, "postFilter");
  ofilt_->postFilter(gv, ov, bv, dv);
  Log::trace() << "ObsFilter<OBS>::postFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
Variables ObsFilter<OBS>::requiredVars() const {
  Log::trace() << "ObsFilter<OBS>::requiredVars" << std::endl;
  return ofilt_->requiredVars();
}

// -----------------------------------------------------------------------------

template <typename OBS>
Variables ObsFilter<OBS>::requiredHdiagnostics() const {
  Log::trace() << "ObsFilter<OBS>::requiredHdiagnostics" << std::endl;
  return ofilt_->requiredHdiagnostics();
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsFilter<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsFilter<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *ofilt_;
  Log::trace() << "ObsFilter<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsFilter<OBS>::checkFilterData(const FilterStage filterStage) {
  Log::trace() << "ObsFilter<OBS>::checkFilterData starting" << std::endl;
  util::Timer timer(classname(), "checkFilterData");
  ofilt_->checkFilterData(filterStage);
  Log::trace() << "ObsFilter<OBS>::checkFilterData done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSFILTER_H_
