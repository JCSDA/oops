/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_HYBRIDGAIN_H_
#define OOPS_RUNS_HYBRIDGAIN_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

class HybridWeightsParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HybridWeightsParameters, Parameters);

 public:
  /// Weight applied to the control member.
  RequiredParameter<double> control{"control", this};

  /// Weight applied to the ensemble mean.
  RequiredParameter<double> ensemble{"ensemble", this};
};

/// Options taken by the HybridGain application.
template <typename MODEL> class HybridGainParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(HybridGainParameters, ApplicationParameters);

 public:
  typedef typename Geometry<MODEL>::Parameters_   GeometryParameters_;
  typedef typename State<MODEL>::Parameters_      StateParameters_;
  typedef typename State<MODEL>::WriteParameters_ StateWriteParameters_;

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Hybrid weights.
  RequiredParameter<HybridWeightsParameters> hybridWeights{"hybrid weights", this};

  /// Hybrid type.
  RequiredParameter<std::string> hybridType{"hybrid type", this};

  /// Control state parameters.
  RequiredParameter<StateParameters_> control{"control", this};

  /// Ensemble mean posterior.
  RequiredParameter<StateParameters_> ensembleMeanPosterior{"ensemble mean posterior", this};

  /// Ensemble mean prior, required if hybrid type is "average increment".
  OptionalParameter<StateParameters_> ensembleMeanPrior{"ensemble mean prior", this};

  /// List of ensemble states.
  RequiredParameter<std::vector<StateParameters_>> ensemble{"ensemble", this};

  /// Output parameter for recentered state
  RequiredParameter<StateWriteParameters_> recenteredOutput{"recentered output", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class HybridGain : public Application {
  typedef Geometry<MODEL>                   Geometry_;
  typedef Increment<MODEL>                  Increment_;
  typedef State<MODEL>                      State_;
  typedef typename State_::Parameters_      StateParameters_;
  typedef typename State_::WriteParameters_ StateWriteParameters_;

  typedef HybridGainParameters<MODEL>  HybridGainParameters_;

 public:
  // -----------------------------------------------------------------------------
  explicit HybridGain(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
  // -----------------------------------------------------------------------------
  virtual ~HybridGain() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    // Deserialize parameters
    HybridGainParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup Geometry
    const Geometry_ resol(params.geometry, this->getComm());

    // Read averaging weights
    const double alphaControl = params.hybridWeights.value().control;
    const double alphaEnsemble = params.hybridWeights.value().ensemble;

    // Read hybrid type
    const std::string hybridType = params.hybridType;

    // Get control state
    const State_ xaControl(resol, params.control);
    Log::test() << "Control: " << std::endl << xaControl << std::endl;
    const Variables vars = xaControl.variables();

    // Get posterior ensemble mean
    const State_ xaEmeanPost(resol, params.ensembleMeanPosterior);
    Log::test() << "Ensemble mean posterior: " << std::endl << xaEmeanPost << std::endl;

    // Compute new center
    State_ xNewCenter(resol, vars, xaControl.validTime());
    if (hybridType == "average analysis") {
      // using average of analysis (following Bonavita)
      // xa_hybrid = a1*xa1 + a2*xa2
      // a1+a2 have to equal to 1

      // Check we don't have an unused option:
      if (params.ensembleMeanPrior.value() != boost::none) {
        throw eckit::BadValue("The HybridGain application expects the option 'ensemble mean prior' "
                              "to be set only when 'hybrid type' is 'average increment', but this "
                              "option was provided with 'average analysis'.");
      }

      ASSERT(alphaControl + alphaEnsemble == 1.0);
      xNewCenter.zero();
      xNewCenter.accumul(alphaControl, xaControl);
      xNewCenter.accumul(alphaEnsemble, xaEmeanPost);
    } else if (hybridType == "average increment") {
      // using average of analysis increments (following Whitaker)
      // xa_hybrid = xf_prior + a1*xinc1 + a2*xinc2
      // Note: a1+a2 no longer need to add to one

      // Check we got the necessary option:
      if (params.ensembleMeanPrior.value() == boost::none) {
        throw eckit::BadValue("The HybridGain application expects the option 'ensemble mean prior' "
                              "to be set when 'hybrid type' is 'average increment', but this "
                              "option was not provided.");
      }

      // Get prior ensemble mean
      const State_ xfEmeanPrior(resol, params.ensembleMeanPrior.value().value());
      Log::test() << "Ensemble mean prior: " << std::endl << xfEmeanPrior << std::endl;
      // compute ensemble mean increment
      Increment_ pertEns(resol, vars, xaControl.validTime());
      pertEns.diff(xaEmeanPost, xfEmeanPrior);
      pertEns *= alphaEnsemble;
      // compute control increment
      Increment_ pertControl(resol, vars, xaControl.validTime());
      pertControl.diff(xaControl, xfEmeanPrior);
      pertControl *= alphaControl;
      // compute hybrid posterior
      xNewCenter = xfEmeanPrior;
      xNewCenter += pertEns;
      xNewCenter += pertControl;
    } else {
      ABORT("Unknown hybrid gain type: " + hybridType);
    }

    // Output new center
    StateWriteParameters_ centeredOutput = params.recenteredOutput;
    centeredOutput.setMember(0);
    xNewCenter.write(centeredOutput);
    Log::test() << "new center : " << xNewCenter << std::endl;

    // Get ensemble parameters
    const std::vector<StateParameters_>& ensParams = params.ensemble;
    const unsigned nens = ensParams.size();

    // Recenter ensemble around new center and save
    for (unsigned jj = 0; jj < nens; ++jj) {
      State_ x(resol, ensParams[jj]);
      Increment_ pert(resol, vars, x.validTime());
      pert.diff(x, xaEmeanPost);
      x = xNewCenter;
      x += pert;

      // Save recentered member
      StateWriteParameters_ recenteredOutput = params.recenteredOutput;
      recenteredOutput.setMember(jj+1);
      x.write(recenteredOutput);
      Log::test() << "Recentered member " << jj << " : " << x << std::endl;
    }

    return 0;
  }
  // -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    HybridGainParameters_ params;
    params.outputSchema(outputPath);
  }
  // -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    HybridGainParameters_ params;
    params.validate(fullConfig);
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::HybridGain<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HYBRIDGAIN_H_
