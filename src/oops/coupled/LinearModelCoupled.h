/*
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/LinearModel.h"
#include "oops/base/LinearModelBase.h"
#include "oops/base/Variables.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Duration.h"
#include "oops/util/Printable.h"

#include "oops/coupled/AuxCoupledModel.h"
#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/IncrementCoupled.h"
#include "oops/coupled/ModelAuxIncrementCoupled.h"
#include "oops/coupled/StateCoupled.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Implementation of a two-model "coupled" linear model. The two linear models run
/// sequentially and are not exchanging any information currently. The two linear
/// models have to use the same time resolution.
template <typename MODEL1, typename MODEL2>
class LinearModelCoupled : public util::Printable {
  typedef GeometryCoupled<MODEL1, MODEL2>           Geometry_;
  typedef IncrementCoupled<MODEL1, MODEL2>          Increment_;
  typedef AuxCoupledModel<MODEL1, MODEL2>           ModelAuxCtl_;
  typedef ModelAuxIncrementCoupled<MODEL1, MODEL2>  ModelAuxInc_;
  typedef StateCoupled<MODEL1, MODEL2>              State_;

 public:
  static const std::string classname() {return "oops::interface::LinearModel";}
  static std::vector<std::string> names() {return {"Coupled"};}

  LinearModelCoupled(const Geometry_ &, const eckit::Configuration &);
  ~LinearModelCoupled() = default;

// Set the linearization trajectory
  void setTrajectory(const State_ &, State_ &, const ModelAuxCtl_ &);

// TL forecast
  void initializeTL(Increment_ &) const;
  void stepTL(Increment_ &, const ModelAuxInc_ &) const;
  void finalizeTL(Increment_ &) const;

// AD forecast
  void initializeAD(Increment_ &) const;
  void stepAD(Increment_ &, ModelAuxInc_ &) const;
  void finalizeAD(Increment_ &) const;

// Information and diagnostics
  const util::Duration & timeResolution() const {return tstep_;}
  const util::Duration & stepTrajectory() const {return trajstep_;}

 private:
  void print(std::ostream &) const override;

  util::Duration tstep_;
  util::Duration trajstep_;
  std::shared_ptr<const Geometry_> geom_;
  std::unique_ptr<LinearModelBase<MODEL1>> tlm1_;
  std::unique_ptr<LinearModelBase<MODEL2>> tlm2_;
  bool parallel_;
};

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
LinearModelCoupled<MODEL1, MODEL2>::LinearModelCoupled(const Geometry_ & geom,
                                                       const eckit::Configuration & config)
  : tstep_(), trajstep_(), geom_(new Geometry_(geom)), tlm1_(), tlm2_(),
    parallel_(geom.isParallel()) {
  Log::trace() << "LinearModelCoupled::LinearModelCoupled starting" << std::endl;
  instantiateLinearModelFactory<MODEL1>();
  instantiateLinearModelFactory<MODEL2>();
  const eckit::LocalConfiguration conf1(config, MODEL1::name());
  const eckit::LocalConfiguration conf2(config, MODEL2::name());
  if (parallel_) {
    if (geom.modelNumber() == 1) {
      tlm1_.reset(LinearModelFactory<MODEL1>::create(geom.geometry1(), conf1));
      tstep_ = tlm1_->timeResolution();
      trajstep_ = tlm1_->stepTrajectory();
    }
    if (geom.modelNumber() == 2) {
      tlm2_.reset(LinearModelFactory<MODEL2>::create(geom.geometry2(), conf2));
      tstep_ = tlm2_->timeResolution();
      trajstep_ = tlm2_->stepTrajectory();
    }
  } else {
    tlm1_.reset(LinearModelFactory<MODEL1>::create(geom.geometry1(), conf1));
    tlm2_.reset(LinearModelFactory<MODEL2>::create(geom.geometry2(), conf2));
    ASSERT(tlm1_->timeResolution() == tlm2_->timeResolution());
    tstep_ = tlm1_->timeResolution();
    ASSERT(tlm1_->stepTrajectory() == tlm2_->stepTrajectory());
    trajstep_ = tlm1_->stepTrajectory();
  }

  Log::trace() << "LinearModelCoupled::LinearModelCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void LinearModelCoupled<MODEL1, MODEL2>::setTrajectory(const State_ & state1,
                                                       State_ & state2,
                                                       const ModelAuxCtl_ & maux) {
  Log::trace() << "LinearModelCoupled::setTrajectory starting" << std::endl;
  if (tlm1_) tlm1_->setTrajectory(state1.state1(), state2.state1(), maux.aux1());
  if (tlm2_) tlm2_->setTrajectory(state1.state2(), state2.state2(), maux.aux2());
  Log::trace() << "LinearModelCoupled::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void LinearModelCoupled<MODEL1, MODEL2>::initializeTL(Increment_ & dx) const {
  Log::trace() << "LinearModelCoupled::initialize starting" << std::endl;
  if (tlm1_) tlm1_->initializeTL(dx.increment1());
  if (tlm2_) tlm2_->initializeTL(dx.increment2());
  Log::trace() << "LinearModelCoupled::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void LinearModelCoupled<MODEL1, MODEL2>::stepTL(Increment_ & dx,
                                                const ModelAuxInc_ & maux) const {
  Log::trace() << "LinearModelCoupled::stepTL starting" << std::endl;
  if (tlm1_) tlm1_->stepTL(dx.increment1(), maux.modelauxincrement1());
  if (tlm2_) tlm2_->stepTL(dx.increment2(), maux.modelauxincrement2());
  Log::trace() << "LinearModelCoupled::stepTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void LinearModelCoupled<MODEL1, MODEL2>::finalizeTL(Increment_ & dx) const {
  Log::trace() << "LinearModelCoupled::finalizeTL starting" << std::endl;
  if (tlm1_) tlm1_->finalizeTL(dx.increment1());
  if (tlm2_) tlm2_->finalizeTL(dx.increment2());
  Log::trace() << "LinearModelCoupled::finalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void LinearModelCoupled<MODEL1, MODEL2>::initializeAD(Increment_ & dx) const {
  Log::trace() << "LinearModelCoupled::initializeAD starting" << std::endl;
  if (tlm1_) tlm1_->initializeAD(dx.increment1());
  if (tlm2_) tlm2_->initializeAD(dx.increment2());
  Log::trace() << "LinearModelCoupled::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void LinearModelCoupled<MODEL1, MODEL2>::stepAD(Increment_ & dx, ModelAuxInc_ & maux) const {
  Log::trace() << "LinearModelCoupled::stepAD starting" << std::endl;
  if (tlm1_) tlm1_->stepAD(dx.increment1(), maux.modelauxincrement1());
  if (tlm2_) tlm2_->stepAD(dx.increment2(), maux.modelauxincrement2());
  Log::trace() << "LinearModelCoupled::stepAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void LinearModelCoupled<MODEL1, MODEL2>::finalizeAD(Increment_ & dx) const {
  Log::trace() << "LinearModelCoupled::finalizeAD starting" << std::endl;
  if (tlm1_) tlm1_->finalizeAD(dx.increment1());
  if (tlm2_) tlm2_->finalizeAD(dx.increment2());
  Log::trace() << "LinearModelCoupled::finalizeAD done" << std::endl;
}


// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void LinearModelCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "LinearModelCoupled::print starting" << std::endl;

  if (parallel_) {
    std::stringstream ss;
    ss.setf(os.flags());
    ss.precision(os.precision());
    if (tlm1_) {
      ss << std::endl << "LinearModelCoupled: " << MODEL1::name() << std::endl;
      ss << *tlm1_ << std::endl;
    }
    if (tlm2_) {
      ss << std::endl << "LinearModelCoupled: " << MODEL2::name() << std::endl;
      ss << *tlm2_ << std::endl;
    }
    util::gatherPrint(os, ss.str(), geom_->getCommPairRanks());
  } else {
    os << std::endl << "LinearModelCoupled: " << MODEL1::name() << std::endl;
    os << *tlm1_ << std::endl;
    os << std::endl << "LinearModelCoupled: " << MODEL2::name() << std::endl;
    os << *tlm2_;
  }

  Log::trace() << "LinearModelCoupled::print done" << std::endl;
}

// -----------------------------------------------------------------------------
}  // namespace oops
