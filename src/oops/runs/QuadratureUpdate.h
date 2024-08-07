/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_QUADRATUREUPDATE_H_
#define OOPS_RUNS_QUADRATUREUPDATE_H_

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/instantiateCostFactory.h"
#include "oops/assimilation/instantiateMinFactory.h"
#include "oops/assimilation/QuadratureSolver.h"
#include "oops/base/Geometry.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/generic/instantiateNormFactory.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/printRunStats.h"

namespace oops {

template <typename MODEL, typename OBS>
class QuadratureUpdateParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(QuadratureUpdateParameters, ApplicationParameters)

 public:
  /// Parameters for ensemble member to be updated.
  RequiredParameter<eckit::LocalConfiguration> ensMemberConfig{"ensemble state", this};

  /// Parameters for ensemble mean.
  RequiredParameter<eckit::LocalConfiguration> ensMeanConfig{"ensemble mean", this};

  /// Parameters for variational assimilation
  RequiredParameter<eckit::LocalConfiguration> varConfig{"variational", this};

  /// Parameters for cost function used in initialization
  RequiredParameter<eckit::LocalConfiguration> finalConfig{"final", this};

  /// Parameters for outputting the analysis increment
  RequiredParameter<eckit::LocalConfiguration> outputConfig{"output", this};

  /// Parameters for quadrature
  RequiredParameter<eckit::LocalConfiguration> quadConfig{"quadrature update", this};

  /// Parameters for cost function used in initialization
  RequiredParameter<eckit::LocalConfiguration> cfConfig{"cost function", this};
};

template <typename MODEL, typename OBS> class QuadratureUpdate : public Application {
  typedef Geometry<MODEL>                   Geometry_;
  typedef ControlVariable<MODEL, OBS>       CtrlVar_;
  typedef ControlIncrement<MODEL, OBS>      CtrlInc_;
  typedef State<MODEL>                      State_;
  typedef State4D<MODEL>                    State4D_;
  typedef ModelAuxControl<MODEL>            ModelAux_;
  typedef ObsAuxControls<OBS>               ObsAux_;
  typedef CostJbTotal<MODEL, OBS>           JbTotal_;
  typedef QuadratureSolver<MODEL, OBS>      QuadSolver_;

  typedef QuadratureUpdateParameters<MODEL, OBS> QuadParams_;

 public:
// -----------------------------------------------------------------------------
  explicit QuadratureUpdate(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateCostFactory<MODEL, OBS>();
    instantiateCovarFactory<MODEL>();
    instantiateMinFactory<MODEL, OBS>();
    instantiateNormFactory<MODEL>();
    instantiateObsErrorFactory<OBS>();
    instantiateObsFilterFactory<OBS>();
    instantiateLinearModelFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~QuadratureUpdate() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    Log::trace() << "QuadratureUpdate: execute start" << std::endl;
    util::printRunStats("QuadratureUpdate start");

//  Deserialize parameters
    QuadParams_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);
    int ens_id = params.ensMemberConfig.value().getInt("member id");

//  Setup cost function
    std::unique_ptr<CostFunction<MODEL, OBS>>
      J(CostFactory<MODEL, OBS>::create(params.cfConfig, this->getComm()));

//  Get auxiliary information for model and obs, needed for constructing control variables
    std::shared_ptr<ObsAux_> oaux   = J->jb().jbObsBias().background();
    std::shared_ptr<ModelAux_> maux = J->jb().jbModBias().background();

//  Setup geometry for reading ensemble data
    const Geometry_ geometry(params.cfConfig.value().getSubConfiguration("geometry"), this->getComm(), mpi::myself());

//  Get the ensemble state and ensemble mean
    std::shared_ptr<State4D_> ens_bg_state = std::make_shared<State4D_>(geometry, params.ensMemberConfig);
    Log::test() << "Background ensemble state: " << *ens_bg_state << std::endl;
    std::shared_ptr<State4D_> mean_bg_state = std::make_shared<State4D_>(geometry, params.ensMeanConfig);
    Log::test() << "Background ensemble mean: " << *mean_bg_state << std::endl;

//  Wrapping the ensemble data in control variables
    CtrlVar_ ens_bg_ctrl(ens_bg_state, maux, oaux);
    CtrlVar_ mean_bg_ctrl(mean_bg_state, maux, oaux);

//  Setup outer loop
    eckit::LocalConfiguration varConf(fullConfig, "variational");
    std::vector<eckit::LocalConfiguration> iterconfs;
    varConf.get("iterations", iterconfs);

/// "Dummy" evaluation of test function, which serves to properly initialize
/// the geometry information. Without this, the program will segfault upon
/// trying to initialize a ControlIncrement.

    PostProcessor<State_> eval_post;
    iterconfs[0].set("linearize", true);
    J->evaluate(mean_bg_ctrl, iterconfs[0], eval_post);
    util::printRunStats("QuadratureUpdate linearize " + std::to_string(0));

//  Taking difference between background and ensemble state to form the background increment.
    CtrlInc_ dx(J->jb());
    dx.diff(ens_bg_ctrl, mean_bg_ctrl);
    
//  Computing the analysis increment via numerical quadrature.
    QuadSolver_ quadsolver(*J);
    quadsolver.solve(dx, params.quadConfig.value());
    
//  Add the analysis increment to the background, forming an analysis ensemble member.
    CtrlVar_ ens_an_ctrl(mean_bg_ctrl);
    J->addIncrement(ens_an_ctrl, dx);

//  Save analysis ensemble member
    eckit::LocalConfiguration outConfig = params.outputConfig.value();
    ens_an_ctrl.states().write(outConfig);  // ACTUAL "WRITE" AT: mpas-jedi/src/mpasjedi/Fields/mpas_field_mod.F90, subroutine write_fields, called in mpas_state_interface_mod.F90 subroutine mpas_state_write_file_c. outConfig passed as the c_conf argument of mpas_state_write_file_c.

    util::printRunStats("QuadratureUpdate end");
    Log::trace() << "QuadratureUpdate: execute done" << std::endl;
    
    return 0;
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    // Note: QuadratureUpdate app doesn't have application level Parameters yet;
    // not validating anything.
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::QuadratureUpdate<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_QUADRATUREUPDATE_H_
