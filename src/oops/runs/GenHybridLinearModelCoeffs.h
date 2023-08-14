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
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"

namespace oops {

/// Options taken by the GenHybridLinearModelCoeffs application.
template <typename MODEL>
class GenHybridLinearModelCoeffsParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(GenHybridLinearModelCoeffsParameters, ApplicationParameters);
  typedef HybridLinearModelParameters<MODEL>       HybridLinearModelParameters_;
  typedef typename Geometry<MODEL>::Parameters_    GeometryParameters_;

 public:
  RequiredParameter<HybridLinearModelParameters_> hybridLinearModel{"hybrid linear model", this};
  RequiredParameter<GeometryParameters_> updateGeometry{"update geometry", this};
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
  : Application(comm) {}

  virtual ~GenHybridLinearModelCoeffs() = default;

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    Parameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    ASSERT(params.hybridLinearModel.value().coeffs.value().output.value() != boost::none);

    const Geometry_ updateGeometry(params.updateGeometry, this->getComm());
    HybridLinearModel_ hybridLinearModel(updateGeometry, params.hybridLinearModel);

    if (!fullConfig.getSubConfiguration("test").empty()) {
      Increment<MODEL> dx(updateGeometry, hybridLinearModel.variables(),
                          params.hybridLinearModel.value().coeffs.value().windowBegin.value());
      dx.ones();
      ModelAuxIncrement<MODEL> mauxinc(updateGeometry, eckit::LocalConfiguration());
      util::DateTime time(params.hybridLinearModel.value().coeffs.value().windowBegin.value());
      while (time < (params.hybridLinearModel.value().coeffs.value().windowBegin.value()
                     + params.hybridLinearModel.value().coeffs.value().windowLength.value())) {
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
