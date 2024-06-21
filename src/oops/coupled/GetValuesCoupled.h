/*
 * (C) Copyright 2022-2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/GetValues.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/StateCoupled.h"
#include "oops/coupled/TraitCoupled.h"
#include "oops/coupled/UtilsCoupled.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"

namespace util {
  class Duration;
}

namespace oops {

/// PreProcessHelper template specialization for coupled traits: delegate to unspecialized template
/// for each model component.
template <typename MODEL1, typename MODEL2>
struct PreProcessHelper<TraitCoupled<MODEL1, MODEL2>> {
  static void preProcessModelData(const State<TraitCoupled<MODEL1, MODEL2>> & state) {
    PreProcessHelper<MODEL1>::preProcessModelData(state.state().state1());
    PreProcessHelper<MODEL2>::preProcessModelData(state.state().state2());
  }
  static void preProcessModelData(const Increment<TraitCoupled<MODEL1, MODEL2>> & increment) {
    PreProcessHelper<MODEL1>::preProcessModelData(increment.increment().increment1());
    PreProcessHelper<MODEL2>::preProcessModelData(increment.increment().increment2());
  }
  static void preProcessModelDataAD(const Increment<TraitCoupled<MODEL1, MODEL2>> & increment) {
    PreProcessHelper<MODEL1>::preProcessModelDataAD(increment.increment().increment1());
    PreProcessHelper<MODEL2>::preProcessModelDataAD(increment.increment().increment2());
  }
};

/// GetValues template specialization for coupled traits
template <typename MODEL1, typename MODEL2, typename OBS>
class GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>:
         private util::ObjectCounter<GetValues<TraitCoupled<MODEL1, MODEL2>, OBS> > {
  typedef Geometry<TraitCoupled<MODEL1, MODEL2>>  Geometry_;
  typedef GeoVaLs<OBS>                            GeoVaLs_;
  typedef Increment<TraitCoupled<MODEL1, MODEL2>> Increment_;
  typedef SampledLocations<OBS>                   SampledLocations_;
  typedef State<TraitCoupled<MODEL1, MODEL2>>     State_;

 public:
  static const std::string classname() {return "oops::GetValues";}

  GetValues(const eckit::Configuration &, const Geometry_ &,
            const util::TimeWindow &,
            const SampledLocations_ &,
            const Variables &, const Variables & varl = Variables());

/// Nonlinear
  void initialize(const util::Duration &);
  void process(const State_ &);
  void finalize();
  void fillGeoVaLs(GeoVaLs_ &);

/// TL
  void initializeTL(const util::Duration &);
  void processTL(const Increment_ &);
  void finalizeTL();
  void fillGeoVaLsTL(GeoVaLs_ &);

/// AD
  void fillGeoVaLsAD(const GeoVaLs_ &);
  void initializeAD();
  void processAD(Increment_ &);
  void finalizeAD(const util::Duration &);

/// Variables that will be required from the State and Increment
  const Variables & linearVariables() const {return linvars_;}
  const Variables & requiredVariables() const {return geovars_;}
  const bool useMethodsTL() const {return false;}

 private:
  const Variables geovars_;   /// Variables needed from both models
  const Variables linvars_;   /// Variables for TL/AD needed from both models
  std::unique_ptr<GetValues<MODEL1, OBS>> getvals1_;
  std::unique_ptr<GetValues<MODEL2, OBS>> getvals2_;
};

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::GetValues(const eckit::Configuration & conf,
                                 const Geometry_ & geom,
                                 const util::TimeWindow & timeWindow,
                                 const SampledLocations_ & locs,
                                 const Variables & vars, const Variables & varl)
  : geovars_(vars), linvars_(varl)
{
  Log::trace() << "GetValues::GetValues start" << std::endl;
  // decide what variables are provided by what model
  std::vector<Variables> splitgeovars = splitVariables(geovars_, geom.geometry().variables());
  std::vector<Variables> splitlinvars = splitVariables(linvars_, geom.geometry().variables());

  if (splitgeovars[0].size() > 0) {
    getvals1_ = std::make_unique<GetValues<MODEL1, OBS>>(conf.getSubConfiguration(MODEL1::name()),
                                 geom.geometry().geometry1(), timeWindow, locs, splitgeovars[0],
                                 splitlinvars[0]);
  }
  if (splitgeovars[1].size() > 0) {
    getvals2_ = std::make_unique<GetValues<MODEL2, OBS>>(conf.getSubConfiguration(MODEL2::name()),
                                 geom.geometry().geometry2(), timeWindow, locs, splitgeovars[1],
                                 splitlinvars[1]);
  }
  Log::trace() << "GetValues::GetValues done" << std::endl;
}

// -----------------------------------------------------------------------------
//  Forward methods (called from nonlinear run)
// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::initialize(const util::Duration & tstep) {
  Log::trace() << "GetValues::initialize start" << std::endl;
  if (getvals1_) getvals1_->initialize(tstep);
  if (getvals2_) getvals2_->initialize(tstep);
  Log::trace() << "GetValues::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::process(const State_ & xx) {
  Log::trace() << "GetValues::process start" << std::endl;
  if (getvals1_) getvals1_->process(xx.state().state1());
  if (getvals2_) getvals2_->process(xx.state().state2());
  Log::trace() << "GetValues::process done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::finalize() {
  Log::trace() << "GetValues::finalize start" << std::endl;
  if (getvals1_) getvals1_->finalize();
  if (getvals2_) getvals2_->finalize();
  Log::trace() << "GetValues::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::fillGeoVaLs(GeoVaLs_ & geovals) {
  Log::trace() << "GetValues::fillGeoVaLs start" << std::endl;

  if (getvals1_) getvals1_->fillGeoVaLs(geovals);
  if (getvals2_) getvals2_->fillGeoVaLs(geovals);

  Log::trace() << "GetValues::fillGeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------
//  TL methods
// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::initializeTL(const util::Duration & tstep) {
  throw eckit::NotImplemented("GetValuesCoupled::initializeTL not implemented", Here());
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::processTL(const Increment_ & dx) {
  throw eckit::NotImplemented("GetValuesCoupled::processTL not implemented", Here());
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::finalizeTL() {
  throw eckit::NotImplemented("GetValuesCoupled::finalizeTL not implemented", Here());
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::fillGeoVaLsTL(GeoVaLs_ & geovals) {
  throw eckit::NotImplemented("GetValuesCoupled::fillGeoVaLsTL not implemented", Here());
}

// -----------------------------------------------------------------------------
//  AD methods
// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::initializeAD() {
  throw eckit::NotImplemented("GetValuesCoupled::initializeAD not implemented", Here());
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::processAD(Increment_ & dx) {
  throw eckit::NotImplemented("GetValuesCoupled::processAD not implemented", Here());
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::finalizeAD(const util::Duration & tstep) {
  throw eckit::NotImplemented("GetValuesCoupled::finalizeAD not implemented", Here());
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::fillGeoVaLsAD(const GeoVaLs_ & geovals) {
  throw eckit::NotImplemented("GetValuesCoupled::fillGeoVaLsAD not implemented", Here());
}

// -----------------------------------------------------------------------------

}  // namespace oops
