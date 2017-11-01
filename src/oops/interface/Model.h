/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_MODEL_H_
#define OOPS_INTERFACE_MODEL_H_

#include <string>

#include <boost/noncopyable.hpp>

#include "util/Logger.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// Encapsulates the nonlinear forecast model

// -----------------------------------------------------------------------------

template <typename MODEL>
class Model : public util::Printable,
              private boost::noncopyable,
              private util::ObjectCounter<Model<MODEL> >  {
  typedef typename MODEL::Model                 Model_;
  typedef Geometry<MODEL>            Geometry_;
  typedef ModelAuxControl<MODEL>      ModelAux_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::Model";}

  Model(const Geometry_ &, const eckit::Configuration &);
  ~Model();

// Run the forecast
  void forecast(State_ &, const ModelAux_ &,
                const util::Duration &, PostProcessor<State_> &) const;

// Information and diagnostics
  const util::Duration & timeResolution() const {return model_->timeResolution();}

 private:
  void initialize(State_ &) const;
  void step(State_ &, const ModelAux_ &) const;
  void finalize(State_ &) const;
  void print(std::ostream &) const;

  boost::scoped_ptr<Model_> model_;
};

// =============================================================================

template<typename MODEL>
Model<MODEL>::Model(const Geometry_ & resol, const eckit::Configuration & conf) : model_() {
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
//  ****** NOTE: The forecast method below is not from MODEL::Model ******
// -----------------------------------------------------------------------------

template<typename MODEL>
void Model<MODEL>::forecast(State_ & xx, const ModelAux_ & maux,
                            const util::Duration & len,
                            PostProcessor<State_> & post) const {
  Log::trace() << "Model<MODEL>::forecast starting" << std::endl;
  util::Timer timer(classname(), "forecast");

  const util::DateTime end(xx.validTime() + len);
  Log::info() << "Model:forecast: forecast starting: " << xx << std::endl;
  this->initialize(xx);
  post.initialize(xx, end, model_->timeResolution());
  while (xx.validTime() < end) {
    this->step(xx, maux);
    post.process(xx);
  }
  post.finalize(xx);
  this->finalize(xx);
  Log::info() << "Model:forecast: forecast finished: " << xx << std::endl;
  ASSERT(xx.validTime() == end);

  Log::trace() << "Model<MODEL>::forecast done" << std::endl;
}

// -----------------------------------------------------------------------------
//  ****** NOTE: The forecast method above is not from MODEL::Model ******
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
void Model<MODEL>::step(State_ & xx, const ModelAux_ & maux) const {
  Log::trace() << "Model<MODEL>::step starting" << std::endl;
  util::Timer timer(classname(), "step");
  model_->step(xx.state(), maux.modelauxcontrol());
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

}  // namespace oops

#endif  // OOPS_INTERFACE_MODEL_H_
