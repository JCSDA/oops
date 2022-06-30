/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_RTPP_H_
#define OOPS_RUNS_RTPP_H_

#include <memory>
#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/StateEnsemble.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// Options taken by the RTPP application.
template <typename MODEL> class RTPPParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(RTPPParameters, ApplicationParameters);

  typedef Geometry<MODEL> Geometry_;
  typedef State<MODEL> State_;

 public:
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef typename State_::WriteParameters_     StateWriteParameters_;
  typedef StateEnsembleParameters<MODEL>        StateEnsembleParameters_;
  RequiredParameter<GeometryParameters_> geometry{
      "geometry", "Geometry parameters", this};
  RequiredParameter<StateEnsembleParameters_> background{
      "background", "Background ensemble states", this};
  RequiredParameter<StateEnsembleParameters_> analysis{
      "analysis", "Analysis ensemble states", this};
  RequiredParameter<float> factor{"factor", "Perturbation factor", this};
  OptionalParameter<Variables> analysisVariables{"analysis variables", this};
  RequiredParameter<StateWriteParameters_> output{
      "output", "analysis mean and ensemble members output", this};
};

/// \brief Application for relaxation to prior perturbation (RTPP) inflation
template <typename MODEL> class RTPP : public Application {
  typedef Geometry<MODEL>                   Geometry_;
  typedef Increment<MODEL>                  Increment_;
  typedef State<MODEL>                      State_;
  typedef StateEnsemble<MODEL>              StateEnsemble_;
  typedef typename State_::WriteParameters_ StateWriteParameters_;

 public:
// -----------------------------------------------------------------------------

  explicit RTPP(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}

// -----------------------------------------------------------------------------

  virtual ~RTPP() = default;

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    RTPPParameters<MODEL> params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup geometry
    const Geometry_ geometry(params.geometry, this->getComm(), oops::mpi::myself());

    const float factor = params.factor.value();

    // Read all ensemble members
    StateEnsemble_ bgens(geometry, params.background);
    StateEnsemble_ anens(geometry, params.analysis);
    const size_t nens = bgens.size();
    ASSERT(nens == anens.size());

    Variables anvars = anens.variables();
    if (params.analysisVariables.value() != boost::none) {
      anvars = params.analysisVariables.value().get();
    }

    // calculate ensemble means
    State_ bg_mean = bgens.mean();
    State_ an_mean = anens.mean();

    Log::test() << "Background member 1:" << bgens[0] << std::endl;
    Log::test() << "Analysis member 1:" << anens[0] << std::endl;

    // update analysis with RTPP
    for (size_t jj = 0; jj < nens; ++jj) {
      // calculate RTPP perturbation
      Increment_ pertTot(geometry, anvars, anens[jj].validTime());
      pertTot.zero();

      Increment_ pert(pertTot);

      pert.diff(bgens[jj], bg_mean);
      pertTot.axpy(factor, pert);

      pert.diff(anens[jj], an_mean);
      pertTot.axpy((1.0 - factor), pert);

      // an_mean contains copies of anens[0] for non-state variables
      // using zero+accumul instead of "=" ensures that only
      // state variables are modified in anens[jj]
      anens[jj].zero();
      anens[jj].accumul(1.0, an_mean);

      // add analysis variable perturbations
      anens[jj] += pertTot;
    }
    Log::test() << "Updated Analysis member 1:" << anens[0] << std::endl;

    // save the analysis mean
    an_mean = anens.mean();   // calculate analysis mean
    Log::test() << "Analysis mean:" << an_mean << std::endl;
    StateWriteParameters_ output = params.output;
    output.setMember(0);
    an_mean.write(output);

    // save the analysis ensemble
    size_t mymember;
    for (size_t jj=0; jj < nens; ++jj) {
      mymember = jj+1;
      output.setMember(mymember);
      anens[jj].write(output);
    }

    return 0;
  }

// -----------------------------------------------------------------------------

  void outputSchema(const std::string & outputPath) const override {
    RTPPParameters<MODEL> params;
    params.outputSchema(outputPath);
  }

// -----------------------------------------------------------------------------

  void validateConfig(const eckit::Configuration & fullConfig) const override {
    RTPPParameters<MODEL> params;
    params.validate(fullConfig);
  }

// -----------------------------------------------------------------------------

 private:
  std::string appname() const override {
    return "oops::RTPP<" + MODEL::name() + ">";
  }

// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_RTPP_H_
