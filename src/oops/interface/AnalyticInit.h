/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_ANALYTICINIT_H_
#define OOPS_INTERFACE_ANALYTICINIT_H_

#include <memory>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Locations.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Timer.h"

namespace oops {

/// \brief Initializes GeoVaLs with analytic formula
template <typename OBS>
class AnalyticInit : private util::ObjectCounter<AnalyticInit<OBS> > {
  typedef typename OBS::AnalyticInit  AnalyticInit_;
  typedef GeoVaLs<OBS>                GeoVaLs_;
  typedef Locations<OBS>              Locations_;

 public:
  static const std::string classname() {return "oops::AnalyticInit";}

  /// constructor (parameters from config)
  explicit AnalyticInit(const eckit::Configuration &);

  /// destructor and copy/move constructors/assignments
  ~AnalyticInit();
  AnalyticInit(const AnalyticInit &) = delete;
  AnalyticInit(AnalyticInit &&) = delete;
  AnalyticInit& operator=(const AnalyticInit &) = delete;
  AnalyticInit& operator=(AnalyticInit &&) = delete;

  /// Fill GeoVaLs with values computed by analytic function
  void fillGeoVaLs(const Locations_ &, GeoVaLs_ &) const;

 private:
  std::unique_ptr<AnalyticInit_> analytic_;
};

// =============================================================================

template<typename OBS>
AnalyticInit<OBS>::AnalyticInit(const eckit::Configuration & conf) {
  Log::trace() << "AnalyticInit<OBS>::AnalyticInit starting" << std::endl;
  util::Timer timer(classname(), "AnalyticInit");
  analytic_.reset(new AnalyticInit_(conf));
  Log::trace() << "AnalyticInit<OBS>::AnalyticInit done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
AnalyticInit<OBS>::~AnalyticInit() {
  Log::trace() << "AnalyticInit<OBS>::~AnalyticInit starting" << std::endl;
  util::Timer timer(classname(), "~AnalyticInit");
  analytic_.reset();
  Log::trace() << "AnalyticInit<OBS>::~AnalyticInit done" << std::endl;
}

// -----------------------------------------------------------------------------
/*! \brief GeoVaLs Analytic Initialization
 *
 * \details **AnalyticInit()** was introduced in May, 2018 (initially as a GeoVaLs
 * constructor) for use with the interpolation test.  The interpolation test
 * requires an initialization of a GeoVaLs object based on the same analytic
 * formulae used for the State initialization (see test::TestStateInterpolation()
 * for further information).  This in turn requires information about the
 * vertical profile in addition to the latitude and longitude positional
 * information in the Locations object.  Currently, this information
 * about the vertical profile is obtained from an existing GeoVaLs object
 * (passed as *gvals*) that represents the output of the State::interpolate()
 * method.  The State.StateGenerate section of the configuration file is
 * also passed to this constructor to provide further information required
 * for the analytic initialization.
 *
 * \date May, 2018: created as a constructor (M. Miesch, JCSDA)
 * \date June, 2018: moved to a method (M. Miesch, JCSDA)
 *
 * \sa test::TestStateInterpolation()
 */
template<typename OBS>
void AnalyticInit<OBS>::fillGeoVaLs(const Locations_ & locs, GeoVaLs_ & gvals) const {
  Log::trace() << "AnalyticInit<OBS>::fillGeoVaLs starting" << std::endl;
  util::Timer timer(classname(), "fillGeoVaLs");
  analytic_->fillGeoVaLs(locs.locations(), gvals.geovals());
  Log::trace() << "AnalyticInit<OBS>::fillGeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_ANALYTICINIT_H_
