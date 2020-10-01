/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2018-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_MODEL_H_
#define OOPS_INTERFACE_MODEL_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/ModelBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// \brief Encapsulates the nonlinear forecast model
/// Note: to see methods that need to be implemented in the forecast model implementation,
/// see ModelBase class.
///
/// Note: implementations of this interface can opt to extract their settings either from
/// a Configuration object or from a subclass of ModelParametersBase.
///
/// In the former case, they should provide a constructor with the following signature:
///
///    Model(const Geometry_ &, const eckit::Configuration &);
///
/// In the latter case, the implementer should first define a subclass of ModelParametersBase
/// holding the settings of the model in question. The implementation of the Model interface
/// should then typedef `Parameters_` to the name of that subclass and provide a constructor with
/// the following signature:
///
///    Model(const Geometry_ &, const Parameters_ &);
// -----------------------------------------------------------------------------

template <typename MODEL>
class Model : public util::Printable,
              private boost::noncopyable,
              private util::ObjectCounter<Model<MODEL> >  {
  typedef ModelBase<MODEL>           ModelBase_;
  typedef Geometry<MODEL>            Geometry_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::Model";}

  Model(const Geometry_ &, const ModelParametersBase &);
  Model(const Geometry_ &, const eckit::Configuration &);
  virtual ~Model();

  /// \brief Run the forecast from state \p xx for \p len time, with \p post postprocessors
  /// Does not need to be implemented in the models
  void forecast(State_ & xx, const ModelAux_ &,
                const util::Duration & len, PostProcessor<State_> & post) const;

  /// \brief Time step for running Model's forecast in oops (frequency with which the
  /// State will be updated)
  const util::Duration & timeResolution() const {return model_->timeResolution();}
  /// \brief Model variables (only used in 4DVar)
  const oops::Variables & variables() const {return model_->variables();}

 private:
  /// \brief Forecast initialization, called before every forecast run
  void initialize(State_ &) const;
  /// \brief Forecast "step", called during forecast run; updates state to the next time
  void step(State_ &, const ModelAux_ &) const;
  /// \brief Forecast finalization; called after each forecast run
  void finalize(State_ &) const;
  /// \brief Print, used in logging
  void print(std::ostream &) const;

  /// \brief Pointer to the Model implementation
  std::unique_ptr<ModelBase_> model_;
};

// =============================================================================

template<typename MODEL>
Model<MODEL>::Model(const Geometry_ & resol, const ModelParametersBase & params)
  : model_()
{
  Log::trace() << "Model<MODEL>::Model starting" << std::endl;
  util::Timer timer(classname(), "Model");
  Log::debug() << "Model config is:" << params << std::endl;
  model_.reset(ModelFactory<MODEL>::create(resol, params));
  Log::trace() << "Model<MODEL>::Model done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Model<MODEL>::Model(const Geometry_ & resol, const eckit::Configuration & conf)
  : Model(resol, validateAndDeserialize<ModelParametersWrapper<MODEL>>(conf).modelParameters)
{}

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
void Model<MODEL>::forecast(State_ & xx, const ModelAux_ & maux,
                            const util::Duration & len,
                            PostProcessor<State_> & post) const {
  Log::trace() << "Model<MODEL>::forecast starting" << std::endl;

  const util::DateTime end(xx.validTime() + len);
  Log::info() << "Model:forecast: forecast starting: " << xx << std::endl;
  this->initialize(xx);
  post.initialize(xx, end, model_->timeResolution());
  post.process(xx);
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
