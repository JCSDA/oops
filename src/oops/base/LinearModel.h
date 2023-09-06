
/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2018-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_LINEARMODEL_H_
#define OOPS_BASE_LINEARMODEL_H_

#include <memory>
#include <string>
#include <utility>

#include <boost/noncopyable.hpp>

#include "oops/base/Geometry.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/generic/LinearModelBase.h"
#include "oops/interface/Increment.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// \brief Abstract linear forecast model used by high level algorithms and applications.
///
/// Note: to see methods that need to be implemented in a generic linear forecast model
/// implementation, see LinearModelBase class in generic/LinearModelBase.h. To see methods that need
/// to be implemented in a MODEL-specific linear forecast model implementation, see
/// interface::LinearModelBase class in interface/LinearModelBase.h.
///
// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearModel : public util::Printable,
                    private boost::noncopyable,
                    private util::ObjectCounter<LinearModel<MODEL> >  {
  typedef LinearModelBase<MODEL>     LinearModelBase_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef ModelAuxControl<MODEL>     ModelAuxCtl_;
  typedef ModelAuxIncrement<MODEL>   ModelAuxInc_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::LinearModel";}

  LinearModel(const Geometry_ &, const eckit::Configuration &);
  virtual ~LinearModel();

  /// \brief Run the linear forecast from increment \p dx for \p len time, with \p post
  /// postprocessors. Does not need to be implemented in the subclasses of LinearModelBase
  void forecastTL(Increment_ & dx, const ModelAuxInc_ &,
                  const util::Duration & len,
                  PostProcessor<Increment_> post = PostProcessor<Increment_>(),
                  PostProcessorTLAD<MODEL> cost = PostProcessorTLAD<MODEL>(),
                  const bool idmodel = false) const;

  /// \brief Run the adjoint linear forecast from increment \p dx for \p len -time, with \p post
  /// postprocessors. Note that the clock ticks backwards. Does not need to be implemented in the
  /// subclasses of LinearModelBase
  void forecastAD(Increment_ & dx, ModelAuxInc_ &,
                  const util::Duration & len,
                  PostProcessor<Increment_> post = PostProcessor<Increment_>(),
                  PostProcessorTLAD<MODEL> cost = PostProcessorTLAD<MODEL>(),
                  const bool idmodel = false) const;

  /// \brief Set the trajectory for the linear model
  void setTrajectory(const State_ &, State_ &, const ModelAuxCtl_ &);

  /// \brief Time step for running LinearModel's forecast in oops (frequency with which the
  /// State will be updated)
  const util::Duration & timeResolution() const {return linearmodel_->timeResolution();}
  /// \brief LinearModel variables (only used in 4DVar)
  const oops::Variables & variables() const {return linearmodel_->variables();}

 private:
  /// \brief Tangent linear forecast initialization, called before every run
  void initializeTL(Increment_ &) const;
  /// \brief Tangent linear forecast "step", called during run; updates increment to the next time
  void stepTL(Increment_ &, const ModelAuxInc_ &) const;
  /// \brief Tangent linear forecast finalization; called after each run
  void finalizeTL(Increment_ &) const;
  /// \brief Adjoint forecast initialization, called before every run
  void initializeAD(Increment_ &) const;
  /// \brief Adjoint forecast "step", called during run; updates increment to the next time
  void stepAD(Increment_ &, ModelAuxInc_ &) const;
  /// \brief Adjoint forecast finalization; called after each run
  void finalizeAD(Increment_ &) const;

  /// \brief Print, used in logging
  void print(std::ostream &) const;

  /// \brief Pointer to the LinearModel implementation
  std::unique_ptr<LinearModelBase_> linearmodel_;
};

// =============================================================================

template<typename MODEL>
LinearModel<MODEL>::LinearModel(const Geometry_ & resol, const eckit::Configuration & config)
  : linearmodel_()
{
  Log::trace() << "LinearModel<MODEL>::LinearModel starting" << std::endl;
  util::Timer timer(classname(), "LinearModel");
  Log::info() << "LinearModel configuration is:" << config << std::endl;
  linearmodel_.reset(LinearModelFactory<MODEL>::create(resol, config));
  Log::trace() << "LinearModel<MODEL>::LinearModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LinearModel<MODEL>::~LinearModel() {
  Log::trace() << "LinearModel<MODEL>::~LinearModel starting" << std::endl;
  util::Timer timer(classname(), "~LinearModel");
  linearmodel_.reset();
  Log::trace() << "LinearModel<MODEL>::~LinearModel done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::forecastTL(Increment_ & dx, const ModelAuxInc_ & maux,
                                    const util::Duration & len,
                                    PostProcessor<Increment_> post,
                                    PostProcessorTLAD<MODEL> cost,
                                    const bool idmodel) const {
  Log::trace() << "LinearModel<MODEL>::forecastTL starting" << std::endl;

  const util::DateTime end(dx.validTime() + len);
  const util::Duration tstep(linearmodel_->timeResolution());
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
      this->stepTL(dx, maux);
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
void LinearModel<MODEL>::forecastAD(Increment_ & dx, ModelAuxInc_ & maux,
                                    const util::Duration & len,
                                    PostProcessor<Increment_> post,
                                    PostProcessorTLAD<MODEL> cost,
                                    const bool idmodel) const {
  Log::trace() << "LinearModel<MODEL>::forecastAD starting" << std::endl;

  const util::DateTime bgn(dx.validTime() - len);
  const util::Duration tstep(linearmodel_->timeResolution());
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
      this->stepAD(dx, maux);
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
void LinearModel<MODEL>::initializeTL(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::initializeTL starting" << std::endl;
  util::Timer timer(classname(), "initializeTL");
  linearmodel_->initializeTL(dx);
  Log::trace() << "LinearModel<MODEL>::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::stepTL(Increment_ & dx, const ModelAuxInc_ & maux) const {
  Log::trace() << "LinearModel<MODEL>::stepTL starting" << std::endl;
  util::Timer timer(classname(), "stepTL");
  linearmodel_->stepTL(dx, maux);
  Log::trace() << "LinearModel<MODEL>::stepTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::finalizeTL(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::finalizeTL starting" << std::endl;
  util::Timer timer(classname(), "finalizeTL");
  linearmodel_->finalizeTL(dx);
  Log::trace() << "LinearModel<MODEL>::finalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::initializeAD(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::initializeAD starting" << std::endl;
  util::Timer timer(classname(), "initializeAD");
  linearmodel_->initializeAD(dx);
  Log::trace() << "LinearModel<MODEL>::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::stepAD(Increment_ & dx, ModelAuxInc_ & maux) const {
  Log::trace() << "LinearModel<MODEL>::stepAD starting" << std::endl;
  util::Timer timer(classname(), "stepAD");
  linearmodel_->stepAD(dx, maux);
  Log::trace() << "LinearModel<MODEL>::stepAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::finalizeAD(Increment_ & dx) const {
  Log::trace() << "LinearModel<MODEL>::finalizeAD starting" << std::endl;
  util::Timer timer(classname(), "finalizeAD");
  linearmodel_->finalizeAD(dx);
  Log::trace() << "LinearModel<MODEL>::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::print(std::ostream & os) const {
  Log::trace() << "LinearModel<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *linearmodel_;
  Log::trace() << "LinearModel<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearModel<MODEL>::setTrajectory(const State_ & xx, State_ & xtraj,
                                       const ModelAuxCtl_ & maux) {
  Log::trace() << "LinearModel<MODEL>::setTrajectory starting" << std::endl;
  util::Timer timer(classname(), "setTrajectory");
  linearmodel_->setTrajectory(xx, xtraj, maux);
  Log::trace() << "LinearModel<MODEL>::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_LINEARMODEL_H_

