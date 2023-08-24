/*
 * (C) Crown copyright 2023, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_RTPS_H_
#define OOPS_RUNS_RTPS_H_

#include <memory>
#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/StateEnsemble.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// Options taken by the RTPS application.
template <typename MODEL> class RTPSParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(RTPSParameters, ApplicationParameters);

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
/// \brief Application for relaxation to prior spread (RTPS) inflation
///
/// \details An application updating the analysis spread using the following equation:
///
/// \f$ xa'_i <- xa'_i(((\alpha * sig_b) + sig_a(1 - \alpha))/sig_a)\f$,
///
/// where sib_b/sig_a refers to background/analysis ensemble standard deviation
/// at each grid point, \alpha is the prescribed factor, and i refers to member i.
///
/// See:
/// Whitaker, J. S., and T. M. Hamill, 2012: Evaluating Methods to Account for
/// System Errors in Ensemble Data Assimilation. Mon. Wea. Rev., 140, 3078–3089,
/// https://doi.org/10.1175/MWR-D-11-00276.1.

template <typename MODEL> class RTPS : public Application {
  typedef Geometry<MODEL>                   Geometry_;
  typedef Increment<MODEL>                  Increment_;
  typedef State<MODEL>                      State_;
  typedef StateEnsemble<MODEL>              StateEnsemble_;
  typedef typename State_::WriteParameters_ StateWriteParameters_;

 public:
// -----------------------------------------------------------------------------

  explicit RTPS(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}

// -----------------------------------------------------------------------------

  virtual ~RTPS() = default;

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    RTPSParameters<MODEL> params;
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

    // calculate ensemble mean
    State_ an_mean = anens.mean();

    // calculate ensemble standard deviations
    Increment_ an_stdDev = anens.variance();
    util::sqrtFieldSet(an_stdDev.fieldSet());
    an_stdDev.synchronizeFields();

    Increment_ bg_stdDev = bgens.variance();
    util::sqrtFieldSet(bg_stdDev.fieldSet());
    bg_stdDev.synchronizeFields();

    // calculate inflation factor
    Increment_ inflation(an_stdDev);
    inflation *= (1.0 - factor);
    inflation.axpy(factor, bg_stdDev);
    util::divideFieldSets(inflation.fieldSet(), an_stdDev.fieldSet());
    inflation.synchronizeFields();

    Log::test() << "Background member 1:" << bgens[0] << std::endl;
    Log::test() << "Analysis member 1:" << anens[0] << std::endl;

    Increment_ pertTot(geometry, anvars, anens[0].validTime());
    // update analysis with RTPS
    for (size_t jj = 0; jj < nens; ++jj) {
      // calculate RTPS perturbation
      pertTot.zero();

      pertTot.diff(anens[jj], an_mean);
      pertTot.schur_product_with(inflation);

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
    RTPSParameters<MODEL> params;
    params.outputSchema(outputPath);
  }

// -----------------------------------------------------------------------------

  void validateConfig(const eckit::Configuration & fullConfig) const override {
    RTPSParameters<MODEL> params;
    params.validate(fullConfig);
  }

// -----------------------------------------------------------------------------

 private:
  std::string appname() const override {
    return "oops::RTPS<" + MODEL::name() + ">";
  }

// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_RTPS_H_
