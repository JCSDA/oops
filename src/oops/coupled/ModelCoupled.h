/*
 * (C) Copyright 2021-2021 UCAR
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

#include "oops/base/Model.h"
#include "oops/base/Variables.h"
#include "oops/generic/ModelBase.h"
#include "oops/interface/ModelBase.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Duration.h"
#include "oops/util/Printable.h"

#include "oops/coupled/AuxCoupledModel.h"
#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/StateCoupled.h"
#include "oops/coupled/TraitCoupled.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Implementation of a two-model "coupled" model. The two models run
/// sequentially and are not exchanging any information currently. The two models
/// have to use the same time resolution.
template <typename MODEL1, typename MODEL2>
class ModelCoupled : public interface::ModelBase<TraitCoupled<MODEL1, MODEL2>> {
  typedef AuxCoupledModel<MODEL1, MODEL2>         AuxCoupledModel_;
  typedef GeometryCoupled<MODEL1, MODEL2>         GeometryCoupled_;
  typedef StateCoupled<MODEL1, MODEL2>            StateCoupled_;

 public:
  ModelCoupled(const GeometryCoupled_ &, const eckit::Configuration &);
  ~ModelCoupled() = default;

  // Run the forecast
  void initialize(StateCoupled_ &) const override;
  void step(StateCoupled_ &, const AuxCoupledModel_ &) const override;
  void finalize(StateCoupled_ &) const override;

  // Information and diagnostics
  const util::Duration & timeResolution() const override {return tstep_;}
  void checkTimes(const StateCoupled_ &) const;

 private:
  void print(std::ostream &) const override;

// Data
  util::Duration tstep_;
  std::shared_ptr<const GeometryCoupled_> geom_;
  std::unique_ptr<ModelBase<MODEL1>> model1_;
  std::unique_ptr<ModelBase<MODEL2>> model2_;
  bool parallel_;
};

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
ModelCoupled<MODEL1, MODEL2>::ModelCoupled(const GeometryCoupled_ & geom,
                                           const eckit::Configuration & config)
  : tstep_(), geom_(new GeometryCoupled_(geom)), model1_(), model2_(),
    parallel_(geom.isParallel()) {
  Log::trace() << "ModelCoupled::ModelCoupled starting" << std::endl;
  Log::debug() << "ModelCoupled::ModelCoupled config " << config << std::endl;
  const eckit::LocalConfiguration conf1(config, MODEL1::name());
  const eckit::LocalConfiguration conf2(config, MODEL2::name());
  if (parallel_) {
    if (geom.modelNumber() == 1) {
      model1_.reset(ModelFactory<MODEL1>::create(geom.geometry1(), conf1));
      tstep_ = model1_->timeResolution();
    }
    if (geom.modelNumber() == 2) {
      model2_.reset(ModelFactory<MODEL2>::create(geom.geometry2(), conf2));
      tstep_ = model2_->timeResolution();
    }
  } else {
    model1_.reset(ModelFactory<MODEL1>::create(geom.geometry1(), conf1));
    model2_.reset(ModelFactory<MODEL2>::create(geom.geometry2(), conf2));
    ASSERT(model1_->timeResolution() == model2_->timeResolution());
    tstep_ = model1_->timeResolution();
  }

  Log::trace() << "ModelCoupled::ModelCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void ModelCoupled<MODEL1, MODEL2>::initialize(StateCoupled_ & xx) const {
  Log::trace() << "ModelCoupled::initialize starting" << std::endl;
  checkTimes(xx);
  if (model1_) model1_->initialize(xx.state1());
  if (model2_) model2_->initialize(xx.state2());
  checkTimes(xx);
  Log::trace() << "ModelCoupled::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void ModelCoupled<MODEL1, MODEL2>::step(StateCoupled_ & xx,
                                        const AuxCoupledModel_ & maux) const {
  Log::trace() << "ModelCoupled::step starting" << std::endl;
  checkTimes(xx);
  if (model1_) model1_->step(xx.state1(), maux.aux1());
  if (model2_) model2_->step(xx.state2(), maux.aux2());
  checkTimes(xx);
  Log::trace() << "ModelCoupled::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void ModelCoupled<MODEL1, MODEL2>::finalize(StateCoupled_ & xx) const {
  Log::trace() << "ModelCoupled::finalize starting" << std::endl;
  checkTimes(xx);
  if (model1_) model1_->finalize(xx.state1());
  if (model2_) model2_->finalize(xx.state2());
  checkTimes(xx);
  Log::trace() << "ModelCoupled::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void ModelCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "ModelCoupled::print starting" << std::endl;

  if (parallel_) {
    std::stringstream ss;
    ss.setf(os.flags());
    ss.precision(os.precision());
    if (model1_) {
      ss << std::endl << "ModelCoupled: " << MODEL1::name() << std::endl;
      ss << *model1_ << std::endl;
    }
    if (model2_) {
      ss << std::endl << "ModelCoupled: " << MODEL2::name() << std::endl;
      ss << *model2_ << std::endl;
    }
    util::gatherPrint(os, ss.str(), geom_->getCommPairRanks());
  } else {
    os << std::endl << "ModelCoupled: " << MODEL1::name() << std::endl;
    os << *model1_ << std::endl;
    os << std::endl << "ModelCoupled: " << MODEL2::name() << std::endl;
    os << *model2_;
  }

  Log::trace() << "ModelCoupled::print done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void ModelCoupled<MODEL1, MODEL2>::checkTimes(const StateCoupled_ & xxs) const {
  if (!parallel_) {
    ASSERT(xxs.state1().validTime() == xxs.state2().validTime());
  } else {
    if (model2_) {
      oops::mpi::send(geom_->getCommPairRanks(), xxs.state2().validTime(), 0, 1234);
    }
    if (model1_) {
      util::DateTime t1(xxs.state1().validTime());
      util::DateTime t2;
      oops::mpi::receive(geom_->getCommPairRanks(), t2, 1, 1234);
      ASSERT(t1 == t2);
    }
  }
}

// -----------------------------------------------------------------------------
}  // namespace oops
