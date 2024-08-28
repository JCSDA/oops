/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSEMBLEGETKFAPPLICATION_H_
#define OOPS_RUNS_ENSEMBLEGETKFAPPLICATION_H_

#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/ForecastParameters.h"
#include "oops/base/Geometry.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/Model.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateSet.h"
#include "oops/base/StateSetSaver.h"
#include "oops/base/StructuredGridPostProcessor.h"

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/mpi/mpi.h"

#include "oops/runs/Application.h"
#include "oops/runs/Forecast.h"

#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the EnsembleForecastApplication application.
template <typename APP>
class EnsembleForecastApplicationParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(EnsembleForecastApplicationParameters, ApplicationParameters)

 public:
  /// Parameters containing a list of YAML files for each ensemble member to be processed.
  RequiredParameter<std::vector<std::string>> files{"files", this};
  RequiredParameter<int> batch{"batch", this};
};

// -----------------------------------------------------------------------------

template <typename APP, typename MODEL>
class EnsembleForecastApplication : public Application {
  typedef EnsembleForecastApplicationParameters<APP> EnsembleForecastApplicationParameters_;
  typedef Geometry<MODEL>              Geometry_;
  typedef Model<MODEL>                 Model_;
  typedef ModelAuxControl<MODEL>       ModelAux_;
  typedef State<MODEL>                 State_;
  typedef StateSet<MODEL>              StateSet_;
  typedef ForecastAppParameters<MODEL> ForecastAppParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit EnsembleForecastApplication(const eckit::mpi::Comm & comm = oops::mpi::world()) :
      Application(comm) {
  }
// -----------------------------------------------------------------------------
  virtual ~EnsembleForecastApplication() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
//  Deserialize parameters
    EnsembleForecastApplicationParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);


//  Get the list of YAML files
    const std::vector<std::string> &files = params.files.value();
    const int batchsize = params.batch.value();
    Log::info() << "params are here " << params << std::endl;
    Log::info() << "EnsembleForecastApplication YAML files:" << files << std::endl;

//  Get the MPI partition
    const int nmembers = files.size();
    const int ntasks = this->getComm().size();
    const int mytask = this->getComm().rank();
    const int tasks_per_member = ntasks / nmembers;
    int mymember = mytask / tasks_per_member + 1;

    Log::info() << "Running " << nmembers << " EnsembleForecastApplication members handled by "
                << ntasks << " total MPI tasks and "
                << tasks_per_member << " MPI tasks per member." << std::endl;

    ASSERT(ntasks%nmembers == 0);

//  Create the communicator for each member, named comm_member_{i}:
    std::string commNameStr = "comm_member_" + std::to_string(mymember);
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commMember = this->getComm().split(mymember, commName);
    const int subrank = commMember.rank();

//  Create the communicator for each face of cubed sphere, named face_member_{i}:
    std::string faceNameStr = "face_member_" + std::to_string(subrank);
    char const *faceName = faceNameStr.c_str();
    eckit::mpi::Comm & faceMember = this->getComm().split(subrank, faceName);
    const int subface = faceMember.rank();

//  Each member uses a different configuration:
    eckit::PathName confPath = files[mymember-1];
    eckit::YAMLConfiguration memberConf(confPath);

    ForecastAppParameters_ fcstparams;
    fcstparams.validate(memberConf);
    fcstparams.deserialize(memberConf);
    const Geometry_ resol(fcstparams.fcstConf.geometry, commMember);

    eckit::LocalConfiguration initialCondition =
          memberConf.getSubConfiguration("initial condition");
    std::vector<int> ens;  // vector of ensemble numbers
    PostProcessor<State_> post;  // Create the post processor where StateSet will be stored
//  Setup times
    Log::info() << "setting up times" << std::endl;
    eckit::LocalConfiguration model = fcstparams.fcstConf.model;
    const util::Duration tstep(model.getString("tstep"));
    eckit::LocalConfiguration ic = fcstparams.fcstConf.initialCondition;
    const util::DateTime bgndate(ic.getString("datetime"));
    const util::Duration fclength = fcstparams.fcstConf.forecastLength;
    const util::DateTime enddate(bgndate + fclength);
    std::vector<util::DateTime> times;

    for (util::DateTime ii=bgndate; ii <= enddate; ii=ii+tstep) {
       times.push_back(ii);
    }
    for (int m = 1; m <=nmembers; m++) { ens.push_back(m); }
    StateSetSaver<MODEL> *saver_ =
          new StateSetSaver<MODEL>(memberConf, resol, times, oops::mpi::myself(), ens, faceMember);
    post.enrollProcessor(saver_);
//  Each member uses a different configuration:
    for (int m = 1; m <=nmembers; m++) {
       if ( m == mymember ) {
         Log::trace() << "running on mymember = " << mymember  << " " << mytask << std::endl;
//         APP ensapp(commMember);
         executeForecast(resol, memberConf, validate, post);
         Log::trace() << "Done with ens execute\n";
       }
       if ( batchsize > 0 ) {  // don't divide by zero
         if (m % batchsize == 0) oops::mpi::world().barrier();
       }
     }
     oops::mpi::world().barrier();
     std::unique_ptr<StateSet_> stateSet = std::move(saver_->getStateSet());

     // calculate the mean and return it in bkg_mean
     StateSet bkg_mean = stateSet->ens_mean();
     Log::trace() << "bkg_mean[0] is " << bkg_mean[0] << std::endl;
     Log::trace() << "bkg_mean[1] is " << bkg_mean[1] << std::endl;
     return 0;
  }
// -----------------------------------------------------------------------------
  void executeForecast(const Geometry_ & resol,
      const eckit::Configuration & fullConfig,
      bool validate, PostProcessor<State_> & post) const {
//  Setup Model
    Log::info() << "Forecast:setting up model" << std::endl;
    const Model_ model(resol, eckit::LocalConfiguration(fullConfig, "model"));

//  Setup initial state
    State_ xx(resol, eckit::LocalConfiguration(fullConfig, "initial condition"));

//  Setup augmented state
    const ModelAux_ moderr(resol, fullConfig.getSubConfiguration("model aux control"));

    const util::Duration fclength(fullConfig.getString("forecast length"));
    const util::DateTime bgndate(xx.validTime());
    const util::DateTime enddate(bgndate + fclength);

    Log::info() << "Forecast:Running forecast from " << bgndate << " to " << enddate << std::endl;
    post.initialize(xx, bgndate, fclength);
//  Run forecast
    Log::info() << "Forecast:running forecast" << std::endl;
    model.forecast(xx, moderr, fclength, post);
    Log::info() << "Forecast:done running forecast" << std::endl;
  }

// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    EnsembleForecastApplicationParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    EnsembleForecastApplicationParameters_ params;
    params.validate(fullConfig);
    // For ensemble applications also need to validate individual yamls
    APP ensapp(oops::mpi::world());
    params.deserialize(fullConfig);
    for (size_t jj = 0; jj < params.files.value().size(); ++jj) {
      const eckit::PathName confPath(params.files.value()[jj]);
      const eckit::YAMLConfiguration memberConf(confPath);
      ensapp.validateConfig(memberConf);
    }
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::EnsembleForecastApplication<>";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ENSEMBLEGETKFAPPLICATION_H_
