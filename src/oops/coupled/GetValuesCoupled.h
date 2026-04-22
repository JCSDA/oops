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
#include "oops/base/Locations.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/StateCoupled.h"
#include "oops/coupled/TraitCoupled.h"
#include "oops/coupled/UtilsCoupled.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"

namespace util {
  class Duration;
}

namespace oops {

/// GetValues template specialization for coupled traits
template <typename MODEL1, typename MODEL2, typename OBS>
class GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>:
         private util::ObjectCounter<GetValues<TraitCoupled<MODEL1, MODEL2>, OBS> > {
  typedef Geometry<TraitCoupled<MODEL1, MODEL2>>  Geometry_;
  typedef GeoVaLs<OBS>                            GeoVaLs_;
  typedef Increment<TraitCoupled<MODEL1, MODEL2>> Increment_;
  typedef Locations<OBS>                          Locations_;
  typedef State<TraitCoupled<MODEL1, MODEL2>>     State_;

 public:
  static const std::string classname() {return "oops::GetValues";}

  GetValues(const eckit::Configuration &, const Geometry_ &,
            const util::TimeWindow &,
            const Locations_ &,
            const Variables &, const Variables & varl = Variables());

  static void preprocess(State_ &);
  static void preprocess(Increment_ &);
  static void preprocessAD(Increment_ &);

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
  bool useMethodsTL() const {return geovalsTL_;}

  /// Continuous DA update
  void updateGetVals(const eckit::Configuration &);

 private:
  const Variables geovars_;   /// Variables needed from both models
  const Variables linvars_;   /// Variables for TL/AD needed from both models
  std::unique_ptr<GetValues<MODEL1, OBS>> getvals1_;
  std::unique_ptr<GetValues<MODEL2, OBS>> getvals2_;
  bool geovalsTL_ = false;
};

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::GetValues(const eckit::Configuration & conf,
                                 const Geometry_ & geom,
                                 const util::TimeWindow & timeWindow,
                                 const Locations_ & locs,
                                 const Variables & vars, const Variables & varl)
  : geovars_(vars), linvars_(varl)
{
  Log::trace() << "GetValuesCoupled::GetValuesCoupled start" << std::endl;
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
  Log::trace() << "GetValuesCoupled::GetValuesCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------
//  Preprocess methods
// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::preprocess(State_ & xx) {
  GetValues<MODEL1, OBS>::preprocess(xx.state().state1());
  GetValues<MODEL2, OBS>::preprocess(xx.state().state2());
}

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::preprocess(Increment_ & dx) {
  GetValues<MODEL1, OBS>::preprocess(dx.increment().increment1());
  GetValues<MODEL2, OBS>::preprocess(dx.increment().increment2());
}

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::preprocessAD(Increment_ & dx) {
  GetValues<MODEL1, OBS>::preprocessAD(dx.increment().increment1());
  GetValues<MODEL2, OBS>::preprocessAD(dx.increment().increment2());
}

// -----------------------------------------------------------------------------
//  Forward methods (called from nonlinear run)
// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::initialize(const util::Duration & tstep) {
  Log::trace() << "GetValuesCoupled::initialize start" << std::endl;
  if (getvals1_) getvals1_->initialize(tstep);
  if (getvals2_) getvals2_->initialize(tstep);
  Log::trace() << "GetValuesCoupled::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::process(const State_ & xx) {
  Log::trace() << "GetValuesCoupled::process start" << std::endl;
  if (getvals1_) getvals1_->process(xx.state().state1());
  if (getvals2_) getvals2_->process(xx.state().state2());
  Log::trace() << "GetValuesCoupled::process done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::finalize() {
  Log::trace() << "GetValuesCoupled::finalize start" << std::endl;
  if (getvals1_) getvals1_->finalize();
  if (getvals2_) getvals2_->finalize();
  Log::trace() << "GetValuesCoupled::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::fillGeoVaLs(GeoVaLs_ & geovals) {
  Log::trace() << "GetValuesCoupled::fillGeoVaLs start" << std::endl;

  if (getvals1_) getvals1_->fillGeoVaLs(geovals);
  if (getvals2_) getvals2_->fillGeoVaLs(geovals);

  Log::trace() << "GetValuesCoupled::fillGeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------
//  TL methods
// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::initializeTL(const util::Duration & tstep) {
  Log::trace() << "GetValuesCoupled::initializeTL start" << std::endl;
  if (getvals1_) getvals1_->initializeTL(tstep);
  if (getvals2_) getvals2_->initializeTL(tstep);
  geovalsTL_ = true;
  Log::trace() << "GetValuesCoupled::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::processTL(const Increment_ & dx) {
  Log::trace() << "GetValuesCoupled::processTL start" << std::endl;
  if (getvals1_) getvals1_->processTL(dx.increment().increment1());
  if (getvals2_) getvals2_->processTL(dx.increment().increment2());
  Log::trace() << "GetValuesCoupled::processTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::finalizeTL() {
  Log::trace() << "GetValuesCoupled::finalizeTL start" << std::endl;
  if (getvals1_) getvals1_->finalizeTL();
  if (getvals2_) getvals2_->finalizeTL();
  Log::trace() << "GetValuesCoupled::finalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::fillGeoVaLsTL(GeoVaLs_ & geovals) {
  Log::trace() << "GetValuesCoupled::fillGeoVaLsTL start" << std::endl;
  if (getvals1_) getvals1_->fillGeoVaLsTL(geovals);
  if (getvals2_) getvals2_->fillGeoVaLsTL(geovals);
  geovalsTL_ = false;
  Log::trace() << "GetValuesCoupled::fillGeoVaLsTL done" << std::endl;
}

// -----------------------------------------------------------------------------
//  AD methods
// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::initializeAD() {
  Log::trace() << "GetValuesCoupled::initializeAD start" << std::endl;
  if (getvals1_) getvals1_->initializeAD();
  if (getvals2_) getvals2_->initializeAD();
  Log::trace() << "GetValuesCoupled::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::processAD(Increment_ & dx) {
  Log::trace() << "GetValuesCoupled::processAD start" << std::endl;
  if (getvals1_) {
    getvals1_->processAD(dx.increment().increment1());
    if (getvals1_->linearVariables().size() > 0) {
      dx.increment().increment1().synchronizeFields();
    }
  }
  if (getvals2_) {
    getvals2_->processAD(dx.increment().increment2());
    if (getvals2_->linearVariables().size() > 0) {
      dx.increment().increment2().synchronizeFields();
    }
  }
  Log::trace() << "GetValuesCoupled::processAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::finalizeAD(const util::Duration & tstep) {
  Log::trace() << "GetValuesCoupled::finalizeAD start" << std::endl;
  if (getvals1_) getvals1_->finalizeAD(tstep);
  if (getvals2_) getvals2_->finalizeAD(tstep);
  Log::trace() << "GetValuesCoupled::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::fillGeoVaLsAD(const GeoVaLs_ & geovals) {
  Log::trace() << "GetValuesCoupled::fillGeoVaLsAD start" << std::endl;
  if (getvals1_) getvals1_->fillGeoVaLsAD(geovals);
  if (getvals2_) getvals2_->fillGeoVaLsAD(geovals);
  Log::trace() << "GetValuesCoupled::fillGeoVaLsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2, typename OBS>
void GetValues<TraitCoupled<MODEL1, MODEL2>, OBS>::updateGetVals(
    const eckit::Configuration & conf) {
  Log::trace() << "GetValuesCoupled::updateGetVals start" << std::endl;
  // passing the same configuration since time window will be the same for both models
  if (getvals1_) getvals1_->updateGetVals(conf);
  if (getvals2_) getvals2_->updateGetVals(conf);
  Log::trace() << "GetValuesCoupled::updateGetVals done" << std::endl;
}

}  // namespace oops
