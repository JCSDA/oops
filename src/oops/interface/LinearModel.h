/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_LINEARMODEL_H_
#define OOPS_INTERFACE_LINEARMODEL_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LinearModelBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/interface/State.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// Encapsulates the linear forecast model.
/*!
 *  This class provides the operations associated with the LinearModel. It wraps 
 *  the actual linear model which can be a model specific one or a generic one 
 *  (identity). The interface for the linear model comprises two levels (LinearModel  
 *  and LinearModelBase) because we want run time polymorphism. 
 *
 *  Note: implementations of this interface can opt to extract their settings either from
 *  a Configuration object or from a subclass of LinearModelParametersBase.
 *
 *  In the former case, they should provide a constructor with the following signature:
 *
 *     LinearModel(const Geometry_ &, const eckit::Configuration &);
 *
 *  In the latter case, the implementer should first define a subclass of LinearModelParametersBase
 *  holding the settings of the model in question. The implementation of the LinearModel interface
 *  should then typedef `Parameters_` to the name of that subclass and provide a constructor with
 *  the following signature:
 *
 *     LinearModel(const Geometry_ &, const Parameters_ &);
 */
// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearModel : public util::Printable,
                    private boost::noncopyable,
                    private util::ObjectCounter<LinearModel<MODEL> >  {
  typedef LinearModelBase<MODEL>     LinearModelBase_;
  typedef Increment<MODEL>           Increment_;
  typedef Geometry<MODEL>            Geometry_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ModelAuxIncrement<MODEL>   ModelAuxIncr_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::LinearModel";}

  LinearModel(const Geometry_ &, const LinearModelParametersBase &);
  LinearModel(const Geometry_ &, const eckit::Configuration &);
  ~LinearModel();

/// Run the tangent linear forecast
  void forecastTL(Increment_ &, const ModelAuxIncr_ &, const util::Duration &,
                  PostProcessor<Increment_> post = PostProcessor<Increment_>(),
                  PostProcessorTLAD<MODEL> cost = PostProcessorTLAD<MODEL>(),
                  const bool idmodel = false) const;

/// Run the adjoint forecast
  void forecastAD(Increment_ &, ModelAuxIncr_ &, const util::Duration &,
                  PostProcessor<Increment_> post = PostProcessor<Increment_>(),
                  PostProcessorTLAD<MODEL> cost = PostProcessorTLAD<MODEL>(),
                  const bool idmodel = false) const;

// Set the linearization trajectory
  void setTrajectory(const State_ &, State_ &, const ModelAux_ &);

// Information and diagnostics
  const util::Duration & timeResolution() const {return tlm_->timeResolution();}
  const oops::Variables & variables() const {return tlm_->variables();}

 protected:
// Run the TL forecast
  void initializeTL(Increment_ &) const;
  void stepTL(Increment_ &, const ModelAuxIncr_ &) const;
  void finalizeTL(Increment_ &) const;

// Run the AD forecast
  void initializeAD(Increment_ &) const;
  void stepAD(Increment_ &, ModelAuxIncr_ &) const;
  void finalizeAD(Increment_ &) const;

 private:
// diagnostics
  void print(std::ostream &) const;

  std::unique_ptr<LinearModelBase_> tlm_;
};

// =============================================================================

template<typename MODEL>
LinearModel<MODEL>::LinearModel(const Geometry_ & resol, const LinearModelParametersBase & params)
  : tlm_()
{
  Log::trace() << "LinearModel<MODEL>::LinearModel starting" << std::endl;
  util::Timer timer(classname(), "LinearModel");
  Log::info() << "LinearModel configuration is:" << params << std::endl;
  tlm_.reset(LinearModelFactory<MODEL>::create(resol, params));
  Log::trace() << "LinearModel<MODEL>::LinearModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LinearModel<MODEL>::LinearModel(const Geometry_ & resol, const eckit::Configuration & conf)
  : LinearModel(resol,
                validateAndDeserialize<LinearModelParametersWrapper<MODEL>>(conf).modelParameters)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
LinearModel<MODEL>::~LinearModel() {
  Log::trace() << "LinearModel<MODEL>::~LinearModel starting" << std::endl;
  util::Timer timer(classname(), "~LinearModel");
  tlm_.reset();
  Log::trace() << "LinearModel<MODEL>::~LinearModel done" << std::endl;
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
/// Run forecast TL and AD
// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::forecastTL(Increment_ & dx, const ModelAuxIncr_ & mctl,
                                    const util::Duration & len,
                                    PostProcessor<Increment_> post,
                                    PostProcessorTLAD<MODEL> cost,
                                    const bool idmodel) const {
  Log::trace() << "LinearModel<MODEL>::forecastTL starting" << std::endl;

  const util::DateTime end(dx.validTime() + len);
  const util::Duration tstep(tlm_->timeResolution());
  Log::info() << "LinearModel<MODEL>::forecastTL: Starting " << dx << std::endl;
  this->initializeTL(dx);
  cost.initializeTL(dx, end, tstep);
  post.initialize(dx, end, tstep);
  cost.processTL(dx);
  post.process(dx);
  if (idmodel) {
    while (dx.validTime() < end) {
      dx.updateTime(tstep);
      cost.processTL(dx);
      post.process(dx);
    }
  } else {
    while (dx.validTime() < end) {
      this->stepTL(dx, mctl);
      cost.processTL(dx);
      post.process(dx);
    }
  }
  cost.finalizeTL(dx);
  post.finalize(dx);
  this->finalizeTL(dx);
  Log::info() << "LinearModel<MODEL>::forecastTL: Finished " << dx << std::endl;
  ASSERT(dx.validTime() == end);

  Log::trace() << "LinearModel<MODEL>::forecastTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::forecastAD(Increment_ & dx, ModelAuxIncr_ & mctl,
                                    const util::Duration & len,
                                    PostProcessor<Increment_> post,
                                    PostProcessorTLAD<MODEL> cost,
                                    const bool idmodel) const {
  Log::trace() << "LinearModel<MODEL>::forecastAD starting" << std::endl;

  const util::DateTime bgn(dx.validTime() - len);
  const util::Duration tstep(tlm_->timeResolution());
  Log::info() << "LinearModel<MODEL>::forecastAD: Starting " << dx << std::endl;
  this->initializeAD(dx);
  post.initialize(dx, bgn, tstep);
  cost.initializeAD(dx, bgn, tstep);
  if (idmodel) {
    while (dx.validTime() > bgn) {
      cost.processAD(dx);
      dx.updateTime(-tstep);
      post.process(dx);
    }
  } else {
    while (dx.validTime() > bgn) {
      cost.processAD(dx);
      this->stepAD(dx, mctl);
      post.process(dx);
    }
  }
  cost.processAD(dx);
  post.process(dx);
  cost.finalizeAD(dx);
  post.finalize(dx);
  this->finalizeAD(dx);
  Log::info() << "LinearModel<MODEL>::forecastAD: Finished " << dx << std::endl;
  ASSERT(dx.validTime() == bgn);

  Log::trace() << "LinearModel<MODEL>::forecastAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::setTrajectory(const State_ & xx, State_ & xlr,
                                       const ModelAux_ & maux) {
  Log::trace() << "LinearModel<MODEL>::setTrajectory starting" << std::endl;
  util::Timer timer(classname(), "setTrajectory");
  tlm_->setTrajectory(xx, xlr, maux);
  Log::trace() << "LinearModel<MODEL>::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::initializeTL(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::initializeTL starting" << std::endl;
  util::Timer timer(classname(), "initializeTL");
  tlm_->initializeTL(dx);
  Log::trace() << "LinearModel<MODEL>::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::stepTL(Increment_ & dx, const ModelAuxIncr_ & merr) const {
  Log::trace() << "LinearModel<MODEL>::stepTL starting" << std::endl;
  util::Timer timer(classname(), "stepTL");
  tlm_->stepTL(dx, merr);
  Log::trace() << "LinearModel<MODEL>::stepTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::finalizeTL(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::finalizeTL starting" << std::endl;
  util::Timer timer(classname(), "finalizeTL");
  tlm_->finalizeTL(dx);
  Log::trace() << "LinearModel<MODEL>::finalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::initializeAD(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::initializeAD starting" << std::endl;
  util::Timer timer(classname(), "initializeAD");
  tlm_->initializeAD(dx);
  Log::trace() << "LinearModel<MODEL>::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::stepAD(Increment_ & dx, ModelAuxIncr_ & merr) const {
  Log::trace() << "LinearModel<MODEL>::stepAD starting" << std::endl;
  util::Timer timer(classname(), "stepAD");
  tlm_->stepAD(dx, merr);
  Log::trace() << "LinearModel<MODEL>::stepAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::finalizeAD(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::finalizeAD starting" << std::endl;
  util::Timer timer(classname(), "finalizeAD");
  tlm_->finalizeAD(dx);
  Log::trace() << "LinearModel<MODEL>::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::print(std::ostream & os) const {
  Log::trace() << "LinearModel<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *tlm_;
  Log::trace() << "LinearModel<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEARMODEL_H_
