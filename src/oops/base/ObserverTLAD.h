/*
 * (C) Copyright 2009-2016 ECMWF.
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
#include "oops/base/GetValueTLAD.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/DateTime.h"

namespace oops {

/// Computes observation equivalent TL and AD to/from increments.

template <typename MODEL, typename OBS>
class ObserverTLAD {
  typedef Geometry<MODEL>              Geometry_;
  typedef GeoVaLs<OBS>                 GeoVaLs_;
  typedef GetValueTLAD<MODEL, OBS>     GetValTLAD_;
  typedef LinearObsOperator<OBS>       LinearObsOperator_;
  typedef Locations<OBS>               Locations_;
  typedef ObsAuxControl<OBS>           ObsAuxCtrl_;
  typedef ObsAuxIncrement<OBS>         ObsAuxIncr_;
  typedef ObsOperator<OBS>             ObsOperator_;
  typedef ObsSpace<OBS>                ObsSpace_;
  typedef ObsVector<OBS>               ObsVector_;

 public:
  ObserverTLAD(const ObsSpace_ &, const eckit::Configuration &);
  ~ObserverTLAD() {}

  std::shared_ptr<GetValTLAD_> initializeTraj(const Geometry_ &, const ObsAuxCtrl_ &);
  void finalizeTraj();

  std::shared_ptr<GetValTLAD_> initializeTL();
  void finalizeTL(const ObsAuxIncr_ &, ObsVector_ &);

  std::shared_ptr<GetValTLAD_> initializeAD(const ObsVector_ &, ObsAuxIncr_ &);
  void finalizeAD();

 private:
  eckit::LocalConfiguration     obsconfig_;
  const ObsSpace_ &             obspace_;    // ObsSpace used in H(x)
  LinearObsOperator_            hoptlad_;    // Linear obs operator
  std::shared_ptr<GetValTLAD_>  getvals_;    // Postproc passed to the model during integration
  std::vector<size_t>  linvars_sizes_;       // Sizes of variables requested from model for
                                             // TL/AD (e.g. number of vertical levels)
  std::unique_ptr<Locations_>   locations_;  // locations
  util::DateTime winbgn_;                    // Begining of assimilation window
  util::DateTime winend_;                    // End of assimilation window
  const ObsAuxCtrl_ *           ybias_;
  bool init_;
};

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
ObserverTLAD<MODEL, OBS>::ObserverTLAD(const ObsSpace_ & obsdb, const eckit::Configuration & conf)
  : obsconfig_(conf), obspace_(obsdb),
    hoptlad_(obspace_,
             validateAndDeserialize<typename LinearObsOperator_::Parameters_>(
               conf.has("linear obs operator") ?
                         eckit::LocalConfiguration(conf, "linear obs operator") :
                         eckit::LocalConfiguration(conf, "obs operator"))),
    getvals_(), locations_(), winbgn_(obsdb.windowStart()), winend_(obsdb.windowEnd()),
    ybias_(nullptr), init_(false)
{
  Log::trace() << "ObserverTLAD::ObserverTLAD" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
std::shared_ptr<GetValueTLAD<MODEL, OBS>>
ObserverTLAD<MODEL, OBS>::initializeTraj(const Geometry_ & geom, const ObsAuxCtrl_ & ybias) {
  Log::trace() << "ObserverTLAD::initializeTraj start" << std::endl;
  ybias_ = &ybias;

//  hop is only needed to get locations and requiredVars
  ObsOperator_ hop(obspace_,
                   validateAndDeserialize<typename ObsOperator_::Parameters_>(
                     eckit::LocalConfiguration(obsconfig_, "obs operator")));
  locations_.reset(new Locations_(hop.locations()));
  linvars_sizes_ = geom.variableSizes(hoptlad_.requiredVars());
  eckit::LocalConfiguration gvconf(obsconfig_.has("linear get values") ?
                         eckit::LocalConfiguration(obsconfig_, "linear get values") :
                           (obsconfig_.has("get values") ?
                            eckit::LocalConfiguration(obsconfig_, "get values") :
                            eckit::LocalConfiguration(obsconfig_, "")));

// Set up variables that will be requested from the model
  Variables geovars;
  geovars += hop.requiredVars();
  geovars += ybias_->requiredVars();

  getvals_.reset(new GetValTLAD_(gvconf, geom, winbgn_, winend_,
                                 *locations_, geovars, hoptlad_.requiredVars()));

  init_ = true;
  Log::trace() << "ObserverTLAD::initializeTraj done" << std::endl;
  return getvals_;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserverTLAD<MODEL, OBS>::finalizeTraj() {
  Log::trace() << "ObserverTLAD::finalizeTraj start" << std::endl;
  ASSERT(init_);

  // GetValues releases GeoVaLs, Observer takes ownership
  std::unique_ptr<GeoVaLs_> geovals = getvals_->finalize();

  /// Set linearization trajectory for H(x)
  hoptlad_.setTrajectory(*geovals, *ybias_);

  init_ = false;
  Log::trace() << "ObserverTLAD::finalizeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
std::shared_ptr<GetValueTLAD<MODEL, OBS>> ObserverTLAD<MODEL, OBS>::initializeTL() {
  Log::trace() << "ObserverTLAD::initializeTL" << std::endl;
  return getvals_;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserverTLAD<MODEL, OBS>::finalizeTL(const ObsAuxIncr_ & ybiastl, ObsVector_ & ydeptl) {
  Log::trace() << "ObserverTLAD::finalizeTL start" << std::endl;

  // GetValues releases GeoVaLs, Observer takes ownership
  std::unique_ptr<GeoVaLs_> geovals = getvals_->finalize();

  // Compute linear H(x)
  hoptlad_.simulateObsTL(*geovals, ydeptl, ybiastl);

  Log::trace() << "ObserverTLAD::finalizeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
std::shared_ptr<GetValueTLAD<MODEL, OBS>>
ObserverTLAD<MODEL, OBS>::initializeAD(const ObsVector_ & ydepad, ObsAuxIncr_ & ybiasad) {
  Log::trace() << "ObserverTLAD::initializeAD start" << std::endl;

  // Compute adjoint of H(x)
  std::unique_ptr<GeoVaLs_> geovals(new GeoVaLs_(*locations_, hoptlad_.requiredVars(),
                                                 linvars_sizes_));
  hoptlad_.simulateObsAD(*geovals, ydepad, ybiasad);

  // GetValues get GeoVaLs and takes ownership
  getvals_->setAD(geovals);

  Log::trace() << "ObserverTLAD::initializeAD done" << std::endl;
  return getvals_;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserverTLAD<MODEL, OBS>::finalizeAD() {
  getvals_->finalizeAD();
  Log::trace() << "ObserverTLAD::finalizeAD done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERTLAD_H_
