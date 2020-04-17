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

#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Increment.h"
#include "oops/interface/LinearGetValues.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

/// Computes observation equivalent TL and AD to/from increments.

template <typename MODEL>
class ObserverTLAD {
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef Increment<MODEL>           Increment_;
  typedef LinearGetValues<MODEL>     LinearGetValues_;
  typedef LinearObsOperator<MODEL>   LinearObsOperator_;
  typedef ObsAuxControl<MODEL>       ObsAuxCtrl_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncr_;
  typedef ObsDiagnostics<MODEL>      ObsDiags_;
  typedef ObsOperator<MODEL>         ObsOperator_;
  typedef ObsSpace<MODEL>            ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;
  typedef State<MODEL>               State_;

 public:
  ObserverTLAD(const eckit::Configuration &,
               const ObsSpace_ &, const ObsAuxCtrl_ &);
  ~ObserverTLAD() {}

  void doInitializeTraj(const State_ &,
                        const util::DateTime &, const util::Duration &, const int &);
  void doProcessingTraj(const State_ &, const util::DateTime &,
                        const util::DateTime &, const int &);
  void doFinalizeTraj(const State_ &);

  void doInitializeTL(const Increment_ &,
                      const util::DateTime &, const util::DateTime &);
  void doProcessingTL(const Increment_ &, const util::DateTime &,
                      const util::DateTime &, const int &);
  void doFinalizeTL(const Increment_ &, ObsVector_ &, const ObsAuxIncr_ &);

  void doFirstAD(Increment_ &, const ObsVector_ &, ObsAuxIncr_ &,
                 const util::DateTime &, const util::DateTime &);
  void doProcessingAD(Increment_ &, const util::DateTime &,
                      const util::DateTime &, const int &);
  void doLastAD(Increment_ &);

 private:
  const ObsSpace_ & obsdb_;
// Obs operator
  ObsOperator_ hop_;
  LinearObsOperator_ hoptlad_;

  const ObsAuxCtrl_ & ybias_;
  Variables geovars_;

  std::vector<std::unique_ptr<LinearGetValues_>> lingetvals_;
  std::shared_ptr<GeoVaLs_> gvals_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObserverTLAD<MODEL>::ObserverTLAD(const eckit::Configuration & config,
                                  const ObsSpace_ & obsdb,
                                  const ObsAuxCtrl_ & ybias)
  : obsdb_(obsdb), hop_(obsdb, eckit::LocalConfiguration(config, "ObsOperator")),
    hoptlad_(obsdb, eckit::LocalConfiguration(config, "LinearObsOperator")),
    ybias_(ybias), geovars_(), lingetvals_(), gvals_()
{
  geovars_ += hop_.variables();
  geovars_ += ybias_.requiredGeoVaLs();
  Log::trace() << "ObserverTLAD::ObserverTLAD" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doInitializeTraj(const State_ & xx,
               const util::DateTime & begin, const util::Duration & winlen, const int & nwin) {
  Log::trace() << "ObserverTLAD::doInitializeTraj start" << std::endl;
  util::DateTime bgn = begin;
  util::DateTime end = bgn + winlen;
  // create one LinearGetValues per subwindow
  for (int iwin = 0; iwin < nwin; ++iwin) {
    lingetvals_.emplace_back(new LinearGetValues_(xx.geometry(),
                                 hop_.locations(bgn, end)));
    bgn = end;
    end = end + winlen;
  }
  gvals_.reset(new GeoVaLs_(hop_.locations(begin, end), geovars_));
  Log::trace() << "ObserverTLAD::doInitializeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doProcessingTraj(const State_ & xx, const util::DateTime & t1,
                                           const util::DateTime & t2, const int & iwin) {
  Log::trace() << "ObserverTLAD::doProcessingTraj start" << std::endl;
// Call nonlinear getValues
  lingetvals_[iwin]->setTrajectory(xx, t1, t2, *gvals_);
  Log::trace() << "ObserverTLAD::doProcessingTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doFinalizeTraj(const State_ & xx) {
  Log::trace() << "ObserverTLAD::doFinalizeTraj start" << std::endl;
  hoptlad_.setTrajectory(*gvals_, ybias_);
  gvals_.reset();
  Log::trace() << "ObserverTLAD::doFinalizeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doInitializeTL(const Increment_ & dx,
                   const util::DateTime & begin, const util::DateTime & end) {
  Log::trace() << "ObserverTLAD::doInitializeTL start" << std::endl;
  gvals_.reset(new GeoVaLs_(hop_.locations(begin, end), hoptlad_.variables()));
  Log::trace() << "ObserverTLAD::doInitializeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doProcessingTL(const Increment_ & dx, const util::DateTime & t1,
                                         const util::DateTime & t2, const int & iwin) {
  Log::trace() << "ObserverTLAD::doProcessingTL start" << std::endl;
// Get increment variables at obs locations
  lingetvals_[iwin]->fillGeoVaLsTL(dx, t1, t2, *gvals_);
  Log::trace() << "ObserverTLAD::doProcessingTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doFinalizeTL(const Increment_ &, ObsVector_ & ydeptl,
                                       const ObsAuxIncr_ & ybiastl) {
  Log::trace() << "ObserverTLAD::doFinalizeTL start" << std::endl;
  hoptlad_.simulateObsTL(*gvals_, ydeptl, ybiastl);
  gvals_.reset();
  Log::trace() << "ObserverTLAD::doFinalizeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doFirstAD(Increment_ & dx, const ObsVector_ & ydepad,
                                    ObsAuxIncr_ & ybiasad,
                                    const util::DateTime & begin,
                                    const util::DateTime & end) {
  Log::trace() << "ObserverTLAD::doFirstAD start" << std::endl;
  gvals_.reset(new GeoVaLs_(hop_.locations(begin, end), hoptlad_.variables()));
  hoptlad_.simulateObsAD(*gvals_, ydepad, ybiasad);
  Log::trace() << "ObserverTLAD::doFirstAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doProcessingAD(Increment_ & dx, const util::DateTime & t1,
                                         const util::DateTime & t2, const int & iwin) {
  Log::trace() << "ObserverTLAD::doProcessingAD start" << std::endl;
// Adjoint of get increment variables at obs locations
  lingetvals_[iwin]->fillGeoVaLsAD(dx, t1, t2, *gvals_);
  Log::trace() << "ObserverTLAD::doProcessingAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObserverTLAD<MODEL>::doLastAD(Increment_ &) {
  Log::trace() << "ObserverTLAD::doLastAD start" << std::endl;
  gvals_.reset();
  Log::trace() << "ObserverTLAD::doLastAD done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERTLAD_H_
