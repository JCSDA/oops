/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "oops/base/ModelBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"

namespace oops {
namespace interface {

// -----------------------------------------------------------------------------

template<typename MODEL>
class Model : public oops::ModelBase<MODEL> {
  typedef typename MODEL::Model         Model_;
  typedef oops::Geometry<MODEL>         Geometry_;
  typedef oops::ModelAuxControl<MODEL>  ModelAux_;
  typedef oops::State<MODEL>            State_;

 public:
  static const std::string classname() {return "oops::interface::Model";}

  Model(const Geometry_ &, const eckit::Configuration &);
  ~Model();

// Forecast
  void initialize(State_ &) const override;
  void step(State_ &, const ModelAux_ &) const override;
  void finalize(State_ &) const override;

// Information and diagnostics
  const util::Duration & timeResolution() const override {return model_->timeResolution();}

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<Model_> model_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
Model<MODEL>::Model(const Geometry_ & resol, const eckit::Configuration & conf)
  : model_()
{
  Log::trace() << "Model<MODEL>::Model starting" << std::endl;
  util::Timer timer(classname(), "Model");
  model_.reset(new Model_(resol.geometry(), conf));
  Log::trace() << "Model<MODEL>::Model done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Model<MODEL>::~Model() {
  Log::trace() << "Model<MODEL>::~Model starting" << std::endl;
  util::Timer timer(classname(), "~Model");
  model_.reset();
  Log::trace() << "Model<MODEL>::~Model done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Model<MODEL>::initialize(State_ & xx) const {
  Log::trace() << "Model<MODEL>::initialize starting" << std::endl;
  util::Timer timer(classname(), "initialize");
  model_->initialize(xx.state());
  Log::trace() << "Model<MODEL>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Model<MODEL>::step(State_ & xx, const ModelAux_ & merr) const {
  Log::trace() << "Model<MODEL>::step starting" << std::endl;
  util::Timer timer(classname(), "step");
  model_->step(xx.state(), merr.modelauxcontrol());
  Log::trace() << "Model<MODEL>::step done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Model<MODEL>::finalize(State_ & xx) const {
  Log::trace() << "Model<MODEL>::finalize starting" << std::endl;
  util::Timer timer(classname(), "finalize");
  model_->finalize(xx.state());
  Log::trace() << "Model<MODEL>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Model<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Model<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *model_;
  Log::trace() << "Model<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace interface
}  // namespace oops
