/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_FINDLOCALOBS_H_
#define OOPS_RUNS_FINDLOCALOBS_H_

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CalcHofX.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class FindLocalObs : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef GeometryIterator<MODEL>    GeometryIterator_;
  typedef Model<MODEL>               Model_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  explicit FindLocalObs(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateObsFilterFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~FindLocalObs() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const eckit::LocalConfiguration windowConf(fullConfig, "Assimilation Window");
    const util::Duration winlen(windowConf.getString("window_length"));
    const util::DateTime winbgn(windowConf.getString("window_begin"));
    const util::DateTime winend(winbgn + winlen);
    Log::info() << "Observation window is:" << windowConf << std::endl;

//  Setup geometry
    const eckit::LocalConfiguration geometryConfig(fullConfig, "Geometry");
    const Geometry_ geometry(geometryConfig, this->getComm());

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "Model");
    const Model_ model(geometry, modelConfig);

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "Initial Condition");
    State_ xx(geometry, model.variables(), initialConfig);
    Log::test() << "Initial state: " << xx << std::endl;

//  Setup forecast outputs
    PostProcessor<State_> post;

    eckit::LocalConfiguration prtConf;
    fullConfig.get("Prints", prtConf);
    post.enrollProcessor(new StateInfo<State_>("fc", prtConf));

//  Setup observations
    eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    ObsSpaces_ obspace(obsconf, this->getComm(), winbgn, winend);

    CalcHofX<MODEL> hofx(obspace, geometry, fullConfig);
    const Observations_ & yobs = hofx.compute(model, xx, post);

//  Save H(x)
    Log::test() << "H(x): " << yobs << std::endl;
    yobs.save("hofx");

//  Iterate over all gridpoints and find local observations (do nothing else for now):
    eckit::LocalConfiguration localconfig(fullConfig, "Localization");
    double dist = localconfig.getDouble("distance");
    int max_nobs = localconfig.getInt("max_nobs");
    for (GeometryIterator_ i = geometry.begin(); i != geometry.end(); ++i) {
       // find all local observations around current gridpoint (*i) with dist and maxnum from config
       ObsSpaces_ localobs(obspace, *i, dist, max_nobs);
       Log::test() << *i << localobs << std::endl;
    }

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::FindLocalObs<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_FINDLOCALOBS_H_
