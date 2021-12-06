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
#include "oops/util/Duration.h"
#include "oops/util/Printable.h"

#include "oops/coupled/AuxCoupledModel.h"
#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/StateCoupled.h"
#include "oops/coupled/TraitCoupled.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Parameters describing a coupled model
template <typename MODEL1, typename MODEL2>
class ModelCoupledParameters : public ModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(ModelCoupledParameters, ModelParametersBase)
  typedef ModelParametersWrapper<MODEL1> Parameters1_;
  typedef ModelParametersWrapper<MODEL2> Parameters2_;
 public:
  RequiredParameter<Parameters1_> model1{MODEL1::name().c_str(), this};
  RequiredParameter<Parameters2_> model2{MODEL2::name().c_str(), this};
};

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
  typedef ModelCoupledParameters<MODEL1, MODEL2>  Parameters_;
  ModelCoupled(const GeometryCoupled_ &, const Parameters_ &);
  ~ModelCoupled() = default;

  // Run the forecast
  void initialize(StateCoupled_ &) const override;
  void step(StateCoupled_ &, const AuxCoupledModel_ &) const override;
  void finalize(StateCoupled_ &) const override;

  // Information and diagnostics
  const util::Duration & timeResolution() const override {return tstep_;}
  const Variables & variables() const override {
    throw eckit::NotImplemented("ModelCoupled::variables", Here());
  }

 private:
  void print(std::ostream &) const override;

// Data
  util::Duration tstep_;
  std::unique_ptr<ModelBase<MODEL1>> model1_;
  std::unique_ptr<ModelBase<MODEL2>> model2_;
};

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
ModelCoupled<MODEL1, MODEL2>::ModelCoupled(const GeometryCoupled_ & geom,
                                           const Parameters_ & params)
  : tstep_(), model1_(), model2_()
{
  Log::trace() << "ModelCoupled::ModelCoupled starting" << std::endl;

  model1_.reset(ModelFactory<MODEL1>::create(geom.geometry1(),
                      params.model1.value().modelParameters));
  model2_.reset(ModelFactory<MODEL2>::create(geom.geometry2(),
                      params.model2.value().modelParameters));

  ASSERT(model1_->timeResolution() == model2_->timeResolution());
  tstep_ = model1_->timeResolution();

  Log::trace() << "ModelCoupled::ModelCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void ModelCoupled<MODEL1, MODEL2>::initialize(StateCoupled_ & xx) const {
  Log::trace() << "ModelCoupled::initialize starting" << std::endl;
  ASSERT(xx.state1().validTime() == xx.state2().validTime());
  model1_->initialize(xx.state1());
  model2_->initialize(xx.state2());
  ASSERT(xx.state1().validTime() == xx.state2().validTime());
  Log::trace() << "ModelCoupled::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void ModelCoupled<MODEL1, MODEL2>::step(StateCoupled_ & xx,
                                             const AuxCoupledModel_ & maux) const {
  Log::trace() << "ModelCoupled::step starting" << std::endl;
  ASSERT(xx.state1().validTime() == xx.state2().validTime());
  model1_->step(xx.state1(), maux.aux1());
  model2_->step(xx.state2(), maux.aux2());
  ASSERT(xx.state1().validTime() == xx.state2().validTime());
  Log::trace() << "ModelCoupled::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void ModelCoupled<MODEL1, MODEL2>::finalize(StateCoupled_ & xx) const {
  Log::trace() << "ModelCoupled::finalize starting" << std::endl;
  ASSERT(xx.state1().validTime() == xx.state2().validTime());
  model1_->finalize(xx.state1());
  model2_->finalize(xx.state2());
  ASSERT(xx.state1().validTime() == xx.state2().validTime());
  Log::trace() << "ModelCoupled::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
void ModelCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "ModelCoupled::print starting" << std::endl;
  os << "ModelCoupled: " << MODEL1::name() << std::endl;
  os << *model1_ << std::endl;
  os << "ModelCoupled: " << MODEL2::name() << std::endl;
  os << *model2_;
  Log::trace() << "ModelCoupled::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

