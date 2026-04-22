/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVERTLAD_H_
#define OOPS_BASE_OBSERVERTLAD_H_

#include <memory>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/GetValues.h"
#include "oops/base/Locations.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/DateTime.h"
#include "oops/util/TimeWindow.h"

namespace oops {

/// Computes observation equivalent TL and AD to/from increments.
template <typename MODEL, typename OBS>
class ObserverTLAD {
  typedef Geometry<MODEL>              Geometry_;
  typedef GeoVaLs<OBS>                 GeoVaLs_;
  typedef GetValues<MODEL, OBS>        GetValues_;
  typedef LinearObsOperator<OBS>       ObsOpTLAD_;
  typedef Locations<OBS>               Locations_;
  typedef ObsAuxControl<OBS>           ObsAuxCtrl_;
  typedef ObsAuxIncrement<OBS>         ObsAuxIncr_;
  typedef ObsOperator<OBS>             ObsOperator_;
  typedef ObsSpace<OBS>                ObsSpace_;
  typedef ObsVector<OBS>               ObsVector_;
  typedef ObsDataVector<OBS, int>      ObsDataInt_;

 public:
  ObserverTLAD(const ObsSpace_ &, const eckit::Configuration &);
  ~ObserverTLAD() {}

  std::shared_ptr<GetValues_> initializeTraj(const Geometry_ &, const ObsAuxCtrl_ &);
  void finalizeTraj(const ObsDataInt_ &);

  void finalizeTL(const ObsAuxIncr_ &, ObsVector_ &);

  void initializeAD(const ObsVector_ &, ObsAuxIncr_ &);
  void finalizeAD() {}

  /// Accessor to linear obs operator
  const ObsOpTLAD_ & linObsOp() {return *hoptlad_;}

 private:
  typedef std::vector<size_t> VariableSizes;

  const ObsSpace_ &           obspace_;
  Variables                   hopVars_;
  VariableSizes               hopVarSizes_;
  Variables                   tladVars_;
  VariableSizes               tladVarSizes_;
  ObsOperator_                hop_;      // obs operator
  std::unique_ptr<ObsOpTLAD_> hoptlad_;  // linear obs operator
  std::shared_ptr<GetValues_> getvals_;
  std::unique_ptr<Locations_> locations_;
  util::TimeWindow            timeWindow_;
  const ObsAuxCtrl_ *         ybias_;
  bool init_;
  eckit::LocalConfiguration gvConf_;
};

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
ObserverTLAD<MODEL, OBS>::ObserverTLAD(const ObsSpace_ & obsdb, const eckit::Configuration & conf)
  : obspace_(obsdb), hopVars_(), hopVarSizes_(),
  tladVars_(), tladVarSizes_(),
  hop_(obspace_, eckit::LocalConfiguration(conf, "obs operator")),
  hoptlad_(), timeWindow_(obsdb.timeWindow()),
  ybias_(nullptr), init_(false),
  gvConf_(conf.getSubConfiguration("get values"))
{
  Log::trace() << "ObserverTLAD::ObserverTLAD" << std::endl;
  if (conf.has("linear obs operator")) {
    hoptlad_ = std::make_unique<ObsOpTLAD_>(obspace_,
                                            eckit::LocalConfiguration(conf, "linear obs operator"));
  } else {
    hoptlad_ = std::make_unique<ObsOpTLAD_>(obspace_,
                                            eckit::LocalConfiguration(conf, "obs operator"));
  }
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
std::shared_ptr<GetValues<MODEL, OBS>>
ObserverTLAD<MODEL, OBS>::initializeTraj(const Geometry_ & geom, const ObsAuxCtrl_ & ybias) {
  Log::trace() << "ObserverTLAD::initializeTraj start" << std::endl;
  ybias_ = &ybias;

  // Get the list of variables to be obtained from the model state
  hopVars_ = hop_.requiredVars();
  hopVars_ += ybias_->requiredVars();
  hopVarSizes_ = geom.variableSizes(hopVars_);

  // Now deal with the variables obtained from the model increment.
  tladVars_ = hoptlad_->requiredVars();
  tladVarSizes_ = geom.variableSizes(tladVars_);

  // Get the observation locations and their discretizations
  locations_ = std::make_unique<Locations_>(hop_.locations());

  // Set up GetValues
  getvals_.reset(new GetValues_(gvConf_, geom, timeWindow_,
                                *locations_, hopVars_, tladVars_));

  init_ = true;
  Log::trace() << "ObserverTLAD::initializeTraj done" << std::endl;
  return getvals_;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserverTLAD<MODEL, OBS>::finalizeTraj(const ObsDataInt_ & qcflags) {
  Log::trace() << "ObserverTLAD::finalizeTraj start" << std::endl;
  ASSERT(init_);

  // Fill geovals
  GeoVaLs_ geovals(*locations_, hopVars_, hopVarSizes_);
  getvals_->fillGeoVaLs(geovals);

  // Compute the reduced representation of the GeoVaLs for which it's been requested
  oops::Variables reducedVars = ybias_->requiredVars();
  hop_.computeReducedVars(reducedVars, geovals);

  /// Set linearization trajectory for H(x)
  hoptlad_->setTrajectory(geovals, *ybias_, qcflags);

  init_ = false;
  Log::trace() << "ObserverTLAD::finalizeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserverTLAD<MODEL, OBS>::finalizeTL(const ObsAuxIncr_ & ybiastl, ObsVector_ & ydeptl) {
  Log::trace() << "ObserverTLAD::finalizeTL start" << std::endl;

  // TODO(wsmigaj): should we allow linear operators to require also *reduced* GeoVaLs?
  // Fill GeoVaLs
  GeoVaLs_ geovals(*locations_, tladVars_, tladVarSizes_);
  getvals_->fillGeoVaLsTL(geovals);

  // Compute linear H(x)
  hoptlad_->simulateObsTL(geovals, ydeptl, ybiastl);

  Log::trace() << "ObserverTLAD::finalizeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserverTLAD<MODEL, OBS>::initializeAD(const ObsVector_ & ydepad, ObsAuxIncr_ & ybiasad) {
  Log::trace() << "ObserverTLAD::initializeAD start" << std::endl;

  GeoVaLs_ geovals(*locations_, tladVars_, tladVarSizes_);

  // (Adjoint of) Compute linear H(x)
  hoptlad_->simulateObsAD(geovals, ydepad, ybiasad);

  // (Adjoint of) Fill geovals
  getvals_->fillGeoVaLsAD(geovals);

  Log::trace() << "ObserverTLAD::initializeAD done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERTLAD_H_
