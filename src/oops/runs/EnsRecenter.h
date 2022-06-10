/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSRECENTER_H_
#define OOPS_RUNS_ENSRECENTER_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// \brief Top-level options taken by the EnsRecenter application.
template <typename MODEL>
class EnsRecenterParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(EnsRecenterParameters, ApplicationParameters)

  typedef Geometry<MODEL> Geometry_;
  typedef State<MODEL> State_;

 public:
  typedef typename Geometry_::Parameters_   GeometryParameters_;
  typedef typename State_::Parameters_      StateParameters_;
  typedef typename State_::WriteParameters_ StateWriteParameters_;

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Central state parameters.
  RequiredParameter<StateParameters_> center{"center", this};

  /// Parameter controlling whether the center should be zeroed out
  Parameter<bool> zeroCenter{"zero center", false, this};

  /// Parameters describing ensemble states to be recentered
  RequiredParameter<std::vector<StateParameters_>> ensemble{"ensemble", this};

  /// Variables to be recentered
  RequiredParameter<oops::Variables> recenterVars{"recenter variables", this};

  /// Parameters for ensemble mean output
  OptionalParameter<StateWriteParameters_> ensmeanOutput{"ensemble mean output", this};

  /// Parameters for recentered ensemble output
  RequiredParameter<StateWriteParameters_> recenteredOutput{"recentered output", this};
};

template <typename MODEL> class EnsRecenter : public Application {
  typedef Geometry<MODEL>   Geometry_;
  typedef Increment<MODEL>  Increment_;
  typedef State<MODEL>      State_;
  typedef typename State_::WriteParameters_ StateWriteParameters_;

 public:
  // -----------------------------------------------------------------------------
  explicit EnsRecenter(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
  // -----------------------------------------------------------------------------
  virtual ~EnsRecenter() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    // Deserialize parameters
    EnsRecenterParameters<MODEL> params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup Geometry
    const Geometry_ resol(params.geometry.value(), this->getComm());

    // Get central state
    State_ x_center(resol, params.center.value());

    // Optionally zero the center
    if (params.zeroCenter) {
      x_center.zero();
    }

    // Get ensemble size
    unsigned nm = params.ensemble.value().size();

    // Compute ensemble mean
    State_ ensmean(x_center);
    ensmean.zero();
    const double rk = 1.0/(static_cast<double>(nm));
    for (unsigned jj = 0; jj < nm; ++jj) {
      State_ x(resol, params.ensemble.value()[jj]);
      ensmean.accumul(rk, x);
      Log::test() << "Original member " << jj << " : " << x << std::endl;
    }
    Log::test() << "Ensemble mean: " << std::endl << ensmean << std::endl;

    // Optionally write the mean out
    if (params.ensmeanOutput.value() != boost::none) {
      ensmean.write(*params.ensmeanOutput.value());
    }

    // Recenter ensemble around central and save
    for (unsigned jj = 0; jj < nm; ++jj) {
      State_ x(resol, params.ensemble.value()[jj]);
      Increment_ pert(resol, params.recenterVars, x.validTime());
      pert.diff(x, ensmean);
      x = x_center;
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
    EnsRecenterParameters<MODEL> params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    EnsRecenterParameters<MODEL> params;
    params.validate(fullConfig);
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::EnsRecenter<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_ENSRECENTER_H_
