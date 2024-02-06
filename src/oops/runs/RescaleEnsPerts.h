/*
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_RESCALEENSPERTS_H_
#define OOPS_RUNS_RESCALEENSPERTS_H_

#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// Options taken by the RescaleEnsPerts application.
template <typename MODEL> class RescaleEnsPertsParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(RescaleEnsPertsParameters, ApplicationParameters);

  typedef Geometry<MODEL> Geometry_;
  typedef Increment<MODEL> Increment_;

 public:
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef typename Increment_::ReadParameters_ IncrementReadParameters_;
  typedef typename Increment_::WriteParameters_ IncrementWriteParameters_;

  RequiredParameter<GeometryParameters_> geometry{
      "geometry", "Geometry parameters", this};
  RequiredParameter<std::vector<IncrementReadParameters_>> sample{
      "sample increments", "Sample of archived analysis increments", this};
  RequiredParameter<std::vector<util::DateTime>> sampleDates{
      "sample dates", "vector of dates corresponding to each increment in the sample", this};
  RequiredParameter<Variables> variables{"variables", this};
  RequiredParameter<IncrementWriteParameters_> output{
      "output", "analysis mean and ensemble members output", this};
  RequiredParameter<util::DateTime> validTime{"valid time", this};
  RequiredParameter<double> factor{"factor", this};
};

// -----------------------------------------------------------------------------
/// Application to carry out the rescaling of ensemble perturbations around the sample mean.
///
/// The application reads in a sample of increments valid at different times.
/// The sample mean is subtracted from each increment and it is then multiplied by
/// the factor given in the configuration.
///
/// This is needed for the method of additive inflation, as described in section 4 of
/// Bowler, N.E., Clayton, A.M., Jardak, M., Lee, E., Lorenc, A.C., Piccolo, C.,
/// Pring, S.R., Wlasak, M.A., Barker, D.M., Inverarity, G.W. and Swinbank, R. (2017),
/// Inflation and localization tests in the development of an ensemble of 4D-ensemble
/// variational assimilations. Q.J.R. Meteorol. Soc., 143: 1280-1302.
/// https://doi.org/10.1002/qj.3004.
///
/// Note that the final step in the paper, of adding the arhcived population mean to the inflated
/// perturbation, is not carried out within this app and is designed to be done by an IAU.
///
template <typename MODEL> class RescaleEnsPerts : public Application {
  typedef Geometry<MODEL>                Geometry_;
  typedef Increment<MODEL>               Increment_;
  typedef IncrementSet<MODEL>            IncrementSet_;

  typedef typename Increment_::ReadParameters_ IncrementReadParameters_;
  typedef typename Increment_::WriteParameters_ IncrementWriteParameters_;

 public:
// -----------------------------------------------------------------------------

  explicit RescaleEnsPerts(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}

// -----------------------------------------------------------------------------

  virtual ~RescaleEnsPerts() {}

// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    Log::trace() << "RescaleEnsPerts: execute start" << std::endl;

    RescaleEnsPertsParameters<MODEL> params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup geometry
    const Geometry_ geometry(params.geometry, this->getComm(), oops::mpi::myself());

    // Setup empty set of increments
    Variables vars = params.variables.value();
    std::vector<util::DateTime> ensTime(1);
    ensTime[0] = params.validTime;
    int rank = params.sample.value().size();
    std::vector<int> ensVector(rank);
    std::iota(ensVector.begin(), ensVector.end(), 0);
    IncrementSet_ increments(geometry, vars, ensTime, oops::mpi::myself(), ensVector);

    // Loop over each increment, temporarily changing the valid times
    // Read in each increment with its corresponding time
    // Update the time back to the valid ensemble time
    for (int jj = 0; jj < rank; ++jj) {
      const util::Duration timeDiff = ensTime[0] - params.sampleDates.value()[jj];
      increments[jj].updateTime(-timeDiff);
      increments[jj].read(params.sample.value()[jj]);
      increments[jj].updateTime(timeDiff);
    }
    Log::test() << "Sample Increments member 1 (time adjusted): " << increments[0] << std::endl;

    IncrementSet_ sampleMean = increments.ens_mean();
    double factor = params.factor.value();
    increments -= sampleMean;
    increments *= factor;

    Log::test() << "Rescaled Perturbation Member 1: " << increments[0] << std::endl;
    for (int jj = 0; jj < rank; ++jj) {
      IncrementWriteParameters_ writeParams = params.output;
      writeParams.setMember(jj);
      increments[jj].write(writeParams);
    }

  return 0;
}
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    RescaleEnsPertsParameters<MODEL> params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    RescaleEnsPertsParameters<MODEL> params;
    params.validate(fullConfig);
  }

// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::RescaleEnsPerts<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_RESCALEENSPERTS_H_
