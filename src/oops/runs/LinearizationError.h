  // (C) Copyright 2023 UCAR.
  // (C) Crown copyright 2023 Met Office.

  // This software is licensed under the terms of the Apache Licence Version 2.0
  // which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#ifndef OOPS_RUNS_LINEARIZATIONERROR_H_
#define OOPS_RUNS_LINEARIZATIONERROR_H_

#include <memory>
#include <string>

#include "oops/base/LinearModel.h"
#include "oops/base/Model.h"
#include "oops/base/StructuredGridWriter.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/generic/instantiateModelFactory.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/printRunStats.h"
#include "oops/util/TimeWindow.h"

namespace oops {

// Requirements for linear model
  // Geometry<MODEL>::Parameters_, "linear model geometry"
  // eckit::LocalConfiguration,    "linear model"
  // eckit::LocalConfiguration,    "model aux increment"

// Requirements for nonlinear model
  // Geometry<MODEL>::Parameters_, "model geometry"
  // eckit::LocalConfiguration,    "model"
  // eckit::LocalConfiguration,    "model aux control"

// Forecast length is the maximum lead time for linearization error computation. Time resolution is
// the difference in lead times (a factor of forecast length; equal to it if just one time desired)
  // util::Duration, "forecastLength"
  // util::Duration, "time resolution"

// Initial nonlinear model states
  // eckit::LocalConfiguration, "x1"
  // eckit::LocalConfiguration, "x2"

// Output configuration for native grid
  // eckit::LocalConfiguration, "output"

// Output configuration for regular structured grid
  // eckit::LocalConfiguration, "output on structured grid"

/// \brief Application for computing and writing fields of linearization error.
///
/// \details Linearization error for a linear model L with respect to a nonlinear model N is defined
/// as:
/// L_[0->T](x_[1,0] - x_[2,0]) - (N_[0->T]x_[1,0] - N_[0->T]x_[2,0])
/// where:
/// L_[0->T] is an operator defining a linear model forecast from time 0 to T,
/// x_[1,0] is some model state at time 0,
/// x_[2,0] is a different model state at time 0,
/// N_[0->T] is an operator defining a nonlinear model forecast from time 0 to T.
///
/// This can be expressed more concisely as:
/// L_[0->T]dx_0 - dx_T
/// where:
/// dx_0 is the difference between the two states at time 0,
/// dx_T is the difference between the two states after being forecasted to time T with N_[0->T].
/// The difference between the two states at time 0 is expected to be suitable for a linear model
/// forecast, e.g. an analysis increment or an ensemble perturbation.
///
/// The equations above assume that both the linear and nonlinear models are on the same grid. This
/// is usually not the case, in which case regridding is required. There are then two options for
/// computing linearization error; on the linear model grid, or on the nonlinear model grid. This
/// application does the latter:
/// P_[l->n]L_[0->T]P_[n->l]dx_0 - dx_T
/// where:
/// P_[l->n] is a regridding from the linear to the nonlinear model grid.
///
/// This application reads in two nonlinear model states, computes the difference between them, and
/// computes linearization error by forecasting the two states with the nonlinear model and the
/// initial difference with the linear model and combining the results using the equations above.
///
/// Linearization error can be computed for a single lead time or for multiple lead times separated
/// by the time resolution specified.
///
/// The computation can be done either at the linear model geometry or the model geometry, and
/// written to file on the native grid, a structured grid (such as lat-lon), or both.
template <typename MODEL>
class LinearizationError : public Application {
  typedef Geometry<MODEL>             Geometry_;
  typedef Increment<MODEL>            Increment_;
  typedef StructuredGridWriter<MODEL>     StructuredGridWriter_;
  typedef LinearModel<MODEL>          LinearModel_;
  typedef Model<MODEL>                Model_;
  typedef ModelAuxControl<MODEL>      ModelAuxControl_;
  typedef ModelAuxIncrement<MODEL>    ModelAuxIncrement_;
  typedef PostProcessorTLAD<MODEL>    PostProcessorTLAD_;
  typedef State<MODEL>                State_;
  typedef TrajectorySaver<MODEL>      TrajectorySaver_;

 public:
  explicit LinearizationError(const eckit::mpi::Comm & comm = oops::mpi::world())
  : Application(comm) {
    instantiateLinearModelFactory<MODEL>();
    instantiateModelFactory<MODEL>();
  }

  virtual ~LinearizationError() = default;

  int execute(const eckit::Configuration & config, bool validate) const override {
    Log::trace() << "LinearizationError::execute() start" << std::endl;
    util::printRunStats("LinearizationError start");

// Set up linear model
    const Geometry_ low(config.getSubConfiguration("linear model geometry"), this->getComm());
    const std::shared_ptr<LinearModel_> linearModel
      = std::make_shared<LinearModel_>(low, config.getSubConfiguration("linear model"));
    const ModelAuxIncrement_ mAuxInc(low, config.getSubConfiguration("model aux increment"));

// Set up model
    const Geometry_ high(config.getSubConfiguration("model geometry"), this->getComm());
    const Model_ model(high, config.getSubConfiguration("model"));
    const ModelAuxControl_ mAuxCtl(high, config.getSubConfiguration("model aux control"));

// Set up StructuredGridWriter
    std::unique_ptr<StructuredGridWriter_> writer;
    if (config.has("output on structured grid")) {
      const eckit::LocalConfiguration structGridConf(config, "output on structured grid");
      writer = std::make_unique<StructuredGridWriter_>(structGridConf, high);
    }

// Set up trajectory saver for use with x1 and empty post processor for x2
    PostProcessor<State_> trajSaver;
    const PostProcessorTLAD_ ppTraj;
    trajSaver.enrollProcessor(new TrajectorySaver_(
      config.getSubConfiguration("linear model"), low, mAuxCtl, linearModel, ppTraj));
    PostProcessor<State_> emptyPp;

// Set up initial states
    State_ x1(high, config.getSubConfiguration("x1"));
    State_ x2(high, config.getSubConfiguration("x2"));
    ASSERT(x1.validTime() == x2.validTime());

// Set up start time, time window, and time resolution at which to compute linearization error
    util::DateTime time(x1.validTime());
    eckit::LocalConfiguration timeWindowConfig;
    timeWindowConfig.set("begin", time.toString());
    timeWindowConfig.set("length", config.getString("forecast length"));
    const util::TimeWindow timeWindow(timeWindowConfig);
    const util::Duration timeResolution(config.getString("time resolution"));
    ASSERT(timeWindow.length() % timeResolution == 0);

// Compute difference between two states, convert to linear model geometry
    Increment_ fdHigh(high, x1.variables(), time);
    fdHigh.diff(x2, x1);
    Increment_ fd(low, fdHigh);
    Log::test() << "fd at " << time << ":" << fd << std::endl;

// Loop over steps of length timeResolution until window is complete
    while (time < timeWindow.end()) {
      time += timeResolution;

      // Forecast states from time to time + timeResolution using model
      model.forecast(x1, mAuxCtl, timeResolution, trajSaver);
      model.forecast(x2, mAuxCtl, timeResolution, emptyPp);

      // Compute difference between two states at time + timeResolution
      fdHigh.updateTime(timeResolution);
      fdHigh.diff(x2, x1);

      // Forecast increment from time to time + timeResolution using linear model
      linearModel->forecastTL(fd, mAuxInc, timeResolution);
      Log::test() << "fd at " << time << ":" << fd << std::endl;

      // Compute linearization error
      Increment_ error(high, fd);
      error -= fdHigh;
      Log::test() << "error at " << time << ":" << error << std::endl;

      // Write to file on structured  grid
      if (writer) {
        writer->interpolateAndWrite(error);
      }
      // Write to file on native grid.
      if (config.has("output")) {
         const eckit::LocalConfiguration outConfig(config, "output");
         error.write(outConfig);
      }
    }

    util::printRunStats("LinearizationError end");
    Log::trace() << "LinearizationError::execute() done" << std::endl;
    return 0;
  }

  void outputSchema(const std::string & outputPath) const override {}

  void validateConfig(const eckit::Configuration & config) const override {}

 private:
  std::string appname() const override {
    return "oops::LinearizationError<" + MODEL::name() + ">";
  }
};

}  // namespace oops

#endif  // OOPS_RUNS_LINEARIZATIONERROR_H_
