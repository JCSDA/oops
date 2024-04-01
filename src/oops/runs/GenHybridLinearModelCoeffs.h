/*
 * (C) Crown copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_GENHYBRIDLINEARMODELCOEFFS_H_
#define OOPS_RUNS_GENHYBRIDLINEARMODELCOEFFS_H_

#include <string>

#include "oops/generic/HybridLinearModel.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/generic/instantiateModelFactory.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"

namespace oops {

/// Options taken by the GenHybridLinearModelCoeffs application.
template <typename MODEL>
class GenHybridLinearModelCoeffsParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(GenHybridLinearModelCoeffsParameters, ApplicationParameters);
 public:
  RequiredParameter<eckit::LocalConfiguration> hybridLinearModel{"hybrid linear model", this};
  RequiredParameter<eckit::LocalConfiguration> updateGeometry{"update geometry", this};
};

/// \brief Application for generating and writing HybridLinearModel coefficients ahead of 4D-Var.
///
/// \details An application that instantiates a HybridLinearModel with the configuration and update
/// geometry specified. Through construction of its HybridLinearModelCoeffs member variable, a set
/// of coefficients is generated and written to disk, for use by a HybridLinearModel in another
/// application.

template <typename MODEL>
class GenHybridLinearModelCoeffs : public Application {
  typedef GenHybridLinearModelCoeffsParameters<MODEL>    Parameters_;
  typedef Geometry<MODEL>                                Geometry_;
  typedef HybridLinearModel<MODEL>                       HybridLinearModel_;

 public:
  explicit GenHybridLinearModelCoeffs(const eckit::mpi::Comm & comm = oops::mpi::world())
  : Application(comm) {
    instantiateLinearModelFactory<MODEL>();
    instantiateModelFactory<MODEL>();
  }

  virtual ~GenHybridLinearModelCoeffs() = default;

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    eckit::LocalConfiguration htlmConf(fullConfig, "hybrid linear model");
    eckit::LocalConfiguration geomConf(fullConfig, "update geometry");
    const Geometry_ updateGeometry(geomConf, this->getComm());
    HybridLinearModel_ hybridLinearModel(updateGeometry, htlmConf);

    if (!fullConfig.getSubConfiguration("test").empty()) {
      const util::TimeWindow timeWindow(htlmConf.getSubConfiguration("coefficients.time window"));
      const Variables vars(fullConfig.getSubConfiguration("test"), "variables");
      Increment<MODEL> dx(updateGeometry, vars, timeWindow.start());
      dx.ones();
      ModelAuxIncrement<MODEL> mauxinc(updateGeometry, eckit::LocalConfiguration());
      util::DateTime time(timeWindow.start());
      while (time < timeWindow.end()) {
        hybridLinearModel.stepTL(dx, mauxinc);
        time += hybridLinearModel.timeResolution();
        Log::test() << "dx at " << time << ": " << dx << std::endl;
      }
    }

    return 0;
  }

  void outputSchema(const std::string & outputPath) const override {
    Parameters_ params;
    params.outputSchema(outputPath);
  }

  void validateConfig(const eckit::Configuration & fullConfig) const override {
    Parameters_ params;
    params.validate(fullConfig);
  }

 private:
  std::string appname() const override {
    return "oops::GenHybridLinearModelCoeffs<" + MODEL::name() + ">";
  }
};

}  // namespace oops

#endif  // OOPS_RUNS_GENHYBRIDLINEARMODELCOEFFS_H_
