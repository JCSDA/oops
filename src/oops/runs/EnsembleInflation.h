/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSEMBLEINFLATION_H_
#define OOPS_RUNS_ENSEMBLEINFLATION_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/Inflation.h"
#include "oops/base/InflationBase.h"
#include "oops/base/instantiateInflationFactory.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

enum class AnalysisType {
  STATE, INCREMENT
};

struct AnalysisTypeParameterTraitsHelper {
  typedef AnalysisType EnumType;
  static constexpr char enumTypeName[] = "AnalysisType";
  static constexpr util::NamedEnumerator<AnalysisType> namedValues[] = {
    { AnalysisType::STATE, "state" },
    { AnalysisType::INCREMENT, "increment" }
  };
};

template <>
struct ParameterTraits<AnalysisType> :
  public EnumParameterTraits<AnalysisTypeParameterTraitsHelper>
{};

constexpr char AnalysisTypeParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<AnalysisType> AnalysisTypeParameterTraitsHelper::namedValues[];

/// Options taken by the EnsembleInflation application.
template <typename MODEL> class EnsembleInflationParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(EnsembleInflationParameters, ApplicationParameters);

  typedef Geometry<MODEL> Geometry_;
  typedef State<MODEL> State_;

 public:
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;

  RequiredParameter<GeometryParameters_> geometry{
      "geometry", "Geometry parameters", this};
  RequiredParameter<eckit::LocalConfiguration> background{
      "background", "Background ensemble states config", this};
  RequiredParameter<eckit::LocalConfiguration> analysis{
      "analysis", "Analysis ensemble config", this};
  OptionalParameter<Variables> analysisVariables{"analysis variables",
      "required input when the analysis type is increment, not required for state", this};
  RequiredParameter<AnalysisType> analysisType{"analysis type", this};
  RequiredParameter<eckit::LocalConfiguration> inflation{"inflation", this};
  RequiredParameter<eckit::LocalConfiguration> output{
      "output", "analysis mean and ensemble members output", this};
};


template <typename MODEL> class EnsembleInflation : public Application {
  typedef Geometry<MODEL>                   Geometry_;
  typedef Increment<MODEL>                  Increment_;
  typedef IncrementSet<MODEL>               IncrementSet_;
  typedef State<MODEL>                      State_;
  typedef StateSet<MODEL>                   StateSet_;


 public:
// -----------------------------------------------------------------------------

  explicit EnsembleInflation(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {
    instantiateInflationFactory<MODEL>();
  }

// -----------------------------------------------------------------------------

  virtual ~EnsembleInflation() = default;

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    EnsembleInflationParameters<MODEL> params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup geometry
    const Geometry_ geometry(params.geometry, this->getComm(), oops::mpi::myself());

    // Read all background ensemble members
    StateSet_ bgens(geometry, params.background);

    // Get inflation subconfigurations
    eckit::LocalConfiguration infConf = params.inflation;
    std::vector<eckit::LocalConfiguration> subconfigs = infConf.getSubConfigurations();

    // Carry out inflation depending on whether analysis is in the form of increments or states
    if (params.analysisType.value() == AnalysisType::STATE) {
      Inflation<MODEL, StateSet_> inflation(params.analysis, geometry, bgens);
      inflation.calculate(subconfigs);
      inflation.save(params.output);
    } else {
      ASSERT(params.analysisVariables.value() != boost::none);
      // Setup analysis variables
      Variables anvars = *params.analysisVariables.value();
      Inflation<MODEL, IncrementSet_> inflation(params.analysis, geometry, bgens, anvars);
      inflation.calculate(subconfigs);
      inflation.save(params.output);
    }
    return 0;
}
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    EnsembleInflationParameters<MODEL> params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    EnsembleInflationParameters<MODEL> params;
    params.validate(fullConfig);
  }

// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::EnsembleInflation<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ENSEMBLEINFLATION_H_
