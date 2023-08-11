/*
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_CONTROLPERT_H_
#define OOPS_RUNS_CONTROLPERT_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/config/YAMLConfiguration.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFctPert.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/IncrementalAssimilation.h"
#include "oops/assimilation/instantiateCostFactory.h"
#include "oops/assimilation/instantiateMinFactory.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/printRunStats.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the ControlPert application.
class ControlPertTemplateParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ControlPertTemplateParameters, ApplicationParameters)

 public:
  /// Template pattern in the \param member template file, to be substituted by the individual
  /// ensemble members' indices left-padded with zeros (subject to \param max number of digits
  /// in member index).
  RequiredParameter<std::string> patternWithPad{"pattern with zero padding", this};

  /// Template pattern in the \param member template file, to be substituted by the individual
  /// ensemble members' indices without zero-padding.
  RequiredParameter<std::string> patternNoPad{"pattern without zero padding", this};

  /// Number of Pert members being run in this executable.
  RequiredParameter<int> nPertMembers{"number of pert members", this};

  /// Member index for the first Pert member that is run in this executable; for example,
  /// if \param nmembers = 4 and \param first pert member index = 17, then the ensemble members
  /// for this run have indices 17, 18, 19 and 20.
  Parameter<int> firstPertMemberIndex{"first pert member index", 1, this};

  /// Maximum number of digits in the ensemble members' indices (the default choice of 3 digits
  /// means that ensemble member indices cannot exceed 999).
  Parameter<int> memberIndexMaxDigits{"max number of digits in member index", 3, this};

  /// Boolean determining whether to run variational DA for the Pert members only.
  Parameter<bool> runPertsOnly{"run pert members only", false, this};
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS> class ControlPert : public Application {
  typedef State<MODEL>                      State_;

 public:
// -----------------------------------------------------------------------------
  explicit ControlPert(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateCostFactory<MODEL, OBS>();
    instantiateCovarFactory<MODEL>();
    instantiateMinFactory<MODEL, OBS>();
    instantiateObsErrorFactory<OBS>();
    instantiateObsFilterFactory<OBS>();
    instantiateLinearModelFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~ControlPert() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    Log::trace() << "ControlPert: execute start" << std::endl;
    util::printRunStats("ControlPert start");

//  Deserialize parameters
    ControlPertTemplateParameters params;
    eckit::LocalConfiguration templateConf(fullConfig, "template");
    if (validate) params.validate(templateConf);
    params.deserialize(templateConf);

//  Retrieve control member configuration (for linearization)
    eckit::LocalConfiguration controlConf(fullConfig, "assimilation");
    const std::string & patternWithPad = params.patternWithPad.value();
    const std::string & patternNoPad = params.patternNoPad.value();
    const int patternLength = params.memberIndexMaxDigits.value();
    util::seekAndReplace(controlConf, patternWithPad, 0, patternLength);
    util::seekAndReplace(controlConf, patternNoPad, std::to_string(0));

//  Get the MPI partition
    const int nPertMembers = params.nPertMembers.value();
    const int nmembers = params.runPertsOnly.value() ? nPertMembers : nPertMembers + 1;
    const int firstPertMember = params.firstPertMemberIndex.value();
    const int lastPertMember = firstPertMember - 1 + nPertMembers;
    const int ntasks = this->getComm().size();
    const int mytask = this->getComm().rank();
    const int tasks_per_member = ntasks / nmembers;
    const int mymember = params.runPertsOnly.value() ?
                             firstPertMember + mytask / tasks_per_member :
                             ((mytask < tasks_per_member) ?
                                   0 : firstPertMember - 1 + mytask / tasks_per_member);
    Log::info() << "Running members "
                << (params.runPertsOnly.value() ? "" : "0 (control member) and ")
                << firstPertMember << " to " << lastPertMember << ", handled by " << ntasks
                << " MPI tasks and " << tasks_per_member << " MPI tasks per member." << std::endl;

    ASSERT(ntasks%nmembers == 0);
    ASSERT(lastPertMember < std::pow(10.0, patternLength));

//  Create the communicator for each member, named comm_member_{i}
    std::string commMemNameStr = "comm_member_" + std::to_string(mymember);
    char const *commMemName = commMemNameStr.c_str();
    eckit::mpi::Comm & commMember = this->getComm().split(mymember, commMemName);

//  Create the communicator for each local area, named comm_area_{i}
    const int myarea = commMember.rank();
    std::string commAreaNameStr = "comm_area_" + std::to_string(myarea);
    char const *commAreaName = commAreaNameStr.c_str();
    eckit::mpi::Comm & commArea = this->getComm().split(myarea, commAreaName);

//  Each member (including control member) uses a different configuration
    eckit::LocalConfiguration memberConf(fullConfig, "assimilation");
    util::seekAndReplace(memberConf, patternWithPad, mymember, patternLength);
    util::seekAndReplace(memberConf, patternNoPad, std::to_string(mymember));

//  Set up cost function for the control member (to be used for the linearization)
    eckit::LocalConfiguration linCostConf(controlConf, "cost function");
    std::unique_ptr<CostFunction<MODEL, OBS>>
      Jlin(CostFactory<MODEL, OBS>::create(linCostConf, commMember));

//  Retrieve cost function type for control member, to determine whether certain variables in
//  the DA for the Pert members need to be time-shifted
    std::string linCostName;
    linCostConf.get("cost type", linCostName);
    const bool shiftTime = !(linCostName == "3D-Var");

//  Check that the control member's cost function is one that uses CostJb3D in the state part
//  of the Jb term, since the Pert members' cost function doesn't currently support other options
//  for CostJbState but it needs to use the same B matrix as the control member's cost function
    ASSERT(linCostName == "3D-Var" || linCostName == "3D-FGAT" || linCostName == "4D-Var");

//  Initialize control member background as linearization state
    ControlVariable<MODEL, OBS> xx(Jlin->jb().getBackground());

//  Set up for the trajectory run
    std::vector<eckit::LocalConfiguration> iterconfsControl;
    controlConf.get("variational.iterations", iterconfsControl);
    PostProcessor<State_> postLin;
    if (iterconfsControl[0].has("prints")) {
      const eckit::LocalConfiguration prtConfig(iterconfsControl[0], "prints");
      postLin.enrollProcessor(new StateInfo<State_>("lintraj", prtConfig));
    }

//  Linearize cost function around linearization state, using configuration for control
//  member's first outer iteration
    iterconfsControl[0].set("iteration", 0);
    iterconfsControl[0].set("linearize", true);
    Jlin->evaluate(xx, iterconfsControl[0], postLin);

//  Retrieve individual members' cost function configurations
    eckit::LocalConfiguration cfConf(memberConf, "cost function");

//  Modify the Pert members' cost function configurations
    if (mymember > 0) {
      // Remove the "cost type" specification
        eckit::LocalConfiguration emptyConf;
        cfConf.set("cost type", emptyConf);

      // Overwrite the "model" sub-configuration (if any) with an empty one
      if (cfConf.has("model")) cfConf.set("model", emptyConf);

      // Remove any prescribed "constraints" terms
      if (cfConf.has("constraints")) {
        std::vector<eckit::LocalConfiguration> emptyConfVector;
        cfConf.set("constraints", emptyConfVector);
      }

      // Set "geometry" to be the first inner-loop resolution of the control member
      eckit::LocalConfiguration geomConf(iterconfsControl[0], "geometry");
      cfConf.set("geometry", geomConf);

      // Turn on observation perturbations (NB: the use of zero-valued observations is guaranteed
      // through the use of the CostJoPert class)
      cfConf.set("observations.obs perturbations", true);

      // Set the nominal "obs perturbations seed" for every ObsSpace to be the member index
      std::vector<eckit::LocalConfiguration> observersConf;
      cfConf.get("observations.observers", observersConf);
      for (unsigned int j = 0; j < observersConf.size(); ++j) {
        observersConf[j].set("obs space.obs perturbations seed", mymember);
      }
      cfConf.set("observations.observers", observersConf);

      // In UFO implementation, no need to remove the "static bc" section under "obs bias", since
      // ObsOperatorPert::simulateObs will not pass any coefficients for static BC predictors to
      // LinearObsOperator::simulateObsTL

      // Update the "cost function" section of the Pert members' configuration accordingly
      memberConf.set("cost function", cfConf);

      Log::info() << "Cost function configuration for Pert member " << mymember << ": "
                  << cfConf << std::endl;
    }

//  Set up cost function for the assimilation
    std::unique_ptr<CostFunction<MODEL, OBS>>
      J((mymember == 0) ? CostFactory<MODEL, OBS>::create(cfConf, commMember) :
        new CostFctPert<MODEL, OBS>(cfConf, commMember, *Jlin, shiftTime));

//  Retrieve length of the assimilation window, to be used if background time-shifting is needed
    std::string linCostWinLength;
    linCostConf.get("window length", linCostWinLength);
    util::Duration winLength(linCostWinLength);

//  Initialize background for all members (for the Pert members, this is the difference between the
//  ensemble member's background and the control member's background, since J->jb().getBackground()
//  has already been modified in the constructor of CostFctPert)
    ControlVariable<MODEL, OBS> xxMember(J->jb().getBackground());

//  Create a copy of xxMember that will not be incremented as the run progresses
    ControlVariable<MODEL, OBS> xxMemberCopy(xxMember);

//  Retrieve individual members' configurations for variational DA
    eckit::LocalConfiguration varConf(memberConf, "variational");

//  Make sure that block minimizers are not used as ControlPert does not support them
    std::string minName;
    varConf.get("minimizer.algorithm", minName);
    ASSERT(minName != "DRPBlockLanczos");

//  Modify the Pert members' variational DA configurations
    if (mymember > 0) {
      // All but the first item in "iterations" are deleted, as ControlPert does not support
      // the running of more than one outer loop in the Pert members
      std::vector<eckit::LocalConfiguration> iterconfsPert;
      varConf.get("iterations", iterconfsPert);
      iterconfsPert.erase(iterconfsPert.begin() + 1, iterconfsPert.end());
      ASSERT(iterconfsPert.size() == 1);

      // Overwrite the "linear model" sub-configuration (if any) with an empty one
      if (iterconfsPert[0].has("linear model")) {
        eckit::LocalConfiguration emptyConf;
        iterconfsPert[0].set("linear model", emptyConf);
      }

      // A switch to let CostFunction::evaluate know that the first-guess gradient of Jb
      // needs to be reset to zero
      iterconfsPert[0].set("control pert", true);

      // Update the "variational" section of the Pert members' configuration accordingly
      varConf.set("iterations", iterconfsPert);
      memberConf.set("variational", varConf);
    }

//  Perform incremental variational assimilation
    int iouter = IncrementalAssimilation<MODEL, OBS>(xxMember, *J, varConf);

//  Retrieve individual members' configurations for the final run of J->evaluate
    eckit::LocalConfiguration finalConfig(memberConf, "final");

//  Set the iteration count (so that the correct QC information can be retrieved from the ODB
//  in ufo::ObsBiasCovariance::linearize), and update the member's configuration accordingly
    finalConfig.set("iteration", iouter);
    memberConf.set("variational", finalConfig);

//  Compute the full analysis increment
    ControlIncrement<MODEL, OBS> dx(J->jb());
    dx.diff(xxMember, xxMemberCopy);

//  Save the full analysis increment if desired
    if (finalConfig.has("increment")) {
      const eckit::LocalConfiguration incConfig(finalConfig, "increment");
      const eckit::LocalConfiguration incGeomConfig(incConfig, "geometry");
      Geometry<MODEL> incGeom(incGeomConfig,
                              dx.state().geometry().getComm(),
                              dx.state().geometry().timeComm());
      Increment<MODEL> dxOutput(incGeom, dx.state());
      const eckit::LocalConfiguration incOutConfig(incConfig, "output");
      dxOutput.write(incOutConfig);
    }

//  Compute final value of the cost function and, if an output configuration is specified,
//  save the control member's analysis trajectory (or state)
    PostProcessor<State_> post;
    if (mymember == 0 && memberConf.has("output")) {
      const eckit::LocalConfiguration outConfig(memberConf, "output");
      post.enrollProcessor(new StateWriter<State_>(outConfig));
    }
    if (finalConfig.has("prints")) {
      const eckit::LocalConfiguration prtConfig(finalConfig, "prints");
      post.enrollProcessor(new StateInfo<State_>("final", prtConfig));
    }
    J->evaluate(xxMember, finalConfig, post);
    Log::info() << "ControlPert: member " << mymember << " incremental assimilation done; "
                << iouter << " iterations." << std::endl;

//  Set up ControlVariable object for the analysis state (valid at the beginning of the
//  assimilation window if the control member's cost function is anything other than 3D-Var;
//  otherwise valid at the middle of the assimilation window)
    ControlVariable<MODEL, OBS> xxMemberAnal(xx);

//  Retrieve the ensemble members' (full-field) background state by adding xxMemberCopy, converted
//  into a ControlIncrement object valid at the same time as xxMember Anal, to xxMemberAnal;
//  then add the Pert increment (time-shifted to the beginning of the window if necessary)
//  to give an intermediate analysis state (which itself has no physical relevance)
    if (mymember > 0) {
      if (shiftTime) xxMemberCopy.state().updateTime(-winLength/2);
      ControlIncrement<MODEL, OBS> xxPertInc(Jlin->jb());
      ControlVariable<MODEL, OBS> xxEmpty(xxMemberCopy, false);
      xxPertInc.diff(xxMemberCopy, xxEmpty);
      J->addIncrement(xxMemberAnal, xxPertInc);
      if (shiftTime) dx.state().updateTime(-winLength/2);
      J->addIncrement(xxMemberAnal, dx);
    }

//  The following computations can only be performed if the control member's DA is run
    if (!params.runPertsOnly.value()) {
//    Add the control member's analysis increment; this returns the analysis state for all
//    members including the control member
      ControlIncrement<MODEL, OBS> dxControl(dx);
      oops::mpi::broadcast(commArea, dxControl, 0);
      J->addIncrement(xxMemberAnal, dxControl);

//    Save the ensemble members' analysis trajectory (or state) if an output configuration
//    is specified
      if (mymember > 0 && memberConf.has("output")) {
        eckit::LocalConfiguration outConfig(memberConf, "output");
        PostProcessor<State_> postMember;
        postMember.enrollProcessor(new StateWriter<State_>(outConfig));
        Jlin->runNL(xxMemberAnal, postMember);  // Hack: by using Jlin instead of J, a model
                                                // trajectory can be computed by calling runNL
                                                // (unless Jlin uses 3D-Var or 4D-Ens-Var)
      }

//    Save the analysis ObsAux for all members including the control member
      xxMemberAnal.obsVar().write(cfConf);
    }

//  Destructing Jlin before J will avoid ObsSpace files associated with J being overwritten,
//  but since Jlin is held in the processor(s) for every ensemble member (including the control
//  member), they need to be destructed one by one to avoid racing condition when saving files
    if (!params.runPertsOnly.value()) {
      if (mymember == 0) Jlin.reset();
      this->getComm().barrier();
    }
    for (int nn = firstPertMember; nn < firstPertMember + nmembers; ++nn) {
      if (mymember == nn) Jlin.reset();
      this->getComm().barrier();
    }

    util::printRunStats("ControlPert end");
    Log::trace() << "ControlPert: execute done" << std::endl;
    return 0;
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    ControlPertTemplateParameters params;
    eckit::LocalConfiguration templateConf(fullConfig, "template");
    params.validate(templateConf);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ControlPert<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_CONTROLPERT_H_
