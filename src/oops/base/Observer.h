/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSERVER_H_
#define OOPS_BASE_OBSERVER_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "oops/base/ObsFilters.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/InterpolatorTraj.h"
#include "oops/interface/LinearObsOperator.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

/// Computes observation equivalent for a single ObsType


// -----------------------------------------------------------------------------

template <typename MODEL>
class Observer : public util::Printable {
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObsDiagnostics<MODEL>      ObsDiags_;
  typedef InterpolatorTraj<MODEL>    InterpolatorTraj_;
  typedef LinearObsOperator<MODEL>   LinearObsOperator_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsAuxControl<MODEL>       ObsAuxCtrl_;
  typedef ObsFilters<MODEL>          ObsFilters_;
  typedef ObsOperator<MODEL>         ObsOperator_;
  typedef ObsVector<MODEL>           ObsVector_;
  typedef State<MODEL>               State_;
  typedef boost::shared_ptr<ObsFilters_> PtrFilters_;

 public:
  Observer(const eckit::Configuration &, const ObsSpace_ &, const ObsAuxCtrl_ &,
           ObsVector_ &,
           const PtrFilters_ filters = PtrFilters_(new ObsFilters_()));
  ~Observer();

  void processTraj(const State_ &, const util::DateTime &, const util::DateTime &,
                   InterpolatorTraj_ &) const;
  void finalizeTraj(const State_ &, LinearObsOperator_ &);

  void doInitialize(const State_ &, const util::DateTime &, const util::DateTime &);
  void doProcessing(const State_ &, const util::DateTime &, const util::DateTime &);
  void doFinalize();

 private:
  void print(std::ostream &) const override;

// Obs operator
  ObsOperator_ hop_;

// Data
  const ObsSpace_ & obsdb_;
  ObsVector_ & yobs_;
  const ObsAuxCtrl_ & ybias_;

  PtrFilters_ filters_;
  Variables geovars_;  // Variables needed from model (through geovals)
  std::shared_ptr<GeoVaLs_> gvals_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Observer<MODEL>::Observer(const eckit::Configuration & conf, const ObsSpace_ & obsdb,
                          const ObsAuxCtrl_ & ybias, ObsVector_ & yobs,
                          const PtrFilters_ filters)
  : hop_(obsdb, eckit::LocalConfiguration(conf, "ObsOperator")),
    obsdb_(obsdb), yobs_(yobs), ybias_(ybias), filters_(filters)
{
  Log::trace() << "Observer::Observer starting" << std::endl;
  geovars_ += hop_.variables();
  geovars_ += ybias_.variables();
  geovars_ += filters_->requiredGeoVaLs();
  Log::trace() << "Observer::Observer done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Observer<MODEL>::~Observer() {
  Log::trace() << "Observer::~Observer starting" << std::endl;
  gvals_.reset();
  Log::trace() << "Observer::~Observer done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Observer<MODEL>::doInitialize(const State_ & xx,
                                   const util::DateTime & begin,
                                   const util::DateTime & end) {
  Log::trace() << "Observer::doInitialize start" << std::endl;
  filters_->preProcess();
  gvals_.reset(new GeoVaLs_(hop_.locations(begin, end), geovars_));
  Log::trace() << "Observer::doInitialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Observer<MODEL>::doProcessing(const State_ & xx,
                                   const util::DateTime & t1,
                                   const util::DateTime & t2) {
  Log::trace() << "Observer::doProcessing start" << std::endl;
// Get state variables at obs locations
  xx.getValues(hop_.locations(t1, t2), geovars_, *gvals_);
  Log::trace() << "Observer::doProcessing done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Observer<MODEL>::processTraj(const State_ & xx,
                                  const util::DateTime & t1,
                                  const util::DateTime & t2,
                                  InterpolatorTraj_ & traj) const {
  Log::trace() << "Observer::processTraj start" << std::endl;
// Get state variables at obs locations and trajectory
  xx.getValues(hop_.locations(t1, t2), geovars_, *gvals_, traj);
  Log::trace() << "Observer::processTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Observer<MODEL>::finalizeTraj(const State_ & xx, LinearObsOperator_ & htlad) {
  Log::trace() << "Observer::finalizeTraj start" << std::endl;
  htlad.setTrajectory(*gvals_, ybias_);
  this->doFinalize();
  Log::trace() << "Observer::finalizeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Observer<MODEL>::doFinalize() {
  Log::trace() << "Observer::doFinalize start" << std::endl;
  filters_->priorFilter(*gvals_);
  ObsDiags_ ydiags(obsdb_, hop_.locations(obsdb_.windowStart(), obsdb_.windowEnd()),
                   filters_->requiredHdiagnostics());
  hop_.simulateObs(*gvals_, yobs_, ybias_, ydiags);
  filters_->postFilter(yobs_, ydiags);
  Log::trace() << "Observer::doFinalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Observer<MODEL>::print(std::ostream &) const {}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVER_H_
