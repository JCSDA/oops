/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_LOCALHOFX_H_
#define OOPS_RUNS_LOCALHOFX_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class LocalHofX : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef GeometryIterator<MODEL>    GeometryIterator_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef ObsAuxControls<MODEL>      ObsAuxCtrls_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  explicit LocalHofX(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateObsFilterFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~LocalHofX() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const eckit::LocalConfiguration windowConf(fullConfig, "Assimilation Window");
    const util::Duration winlen(windowConf.getString("Length"));
    const util::DateTime winbgn(windowConf.getString("Begin"));
    const util::DateTime winend(winbgn + winlen);
    Log::info() << "Observation window is:" << windowConf << std::endl;

//  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "Geometry");
    const Geometry_ resol(resolConfig, this->getComm());

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "Model");
    const Model_ model(resol, modelConfig);

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "Initial Condition");
    Log::info() << "Initial configuration is:" << initialConfig << std::endl;
    State_ xx(resol, model.variables(), initialConfig);
    Log::test() << "Initial state: " << xx << std::endl;

//  Setup augmented state
    ModelAux_ moderr(resol, initialConfig);

//  Setup forecast outputs
    PostProcessor<State_> post;

    eckit::LocalConfiguration prtConf;
    fullConfig.get("Prints", prtConf);
    post.enrollProcessor(new StateInfo<State_>("fc", prtConf));

//  Setup observations
    eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    Log::debug() << "Observations configuration is:" << obsconf << std::endl;
    ObsSpaces_ obsdb(obsconf, this->getComm(), winbgn, winend);

//  Read localization parameters
    eckit::LocalConfiguration localconfig(fullConfig, "Localization");
    double dist = localconfig.getDouble("distance");
    int max_nobs = localconfig.getInt("max_nobs");

//  Get points for finding local obs
    std::vector<eckit::LocalConfiguration> centerconf;
    fullConfig.get("GeoLocations", centerconf);
    std::vector<eckit::geometry::Point2> centers;
    std::vector<boost::shared_ptr<ObsSpaces_>> localobs;
    std::vector<boost::shared_ptr<ObsAuxCtrls_>> localobias;
    std::vector<boost::shared_ptr<Observers<MODEL, State_> >> pobs;
    for (std::size_t jj = 0; jj < centerconf.size(); ++jj) {
       double lon = centerconf[jj].getDouble("lon");
       double lat = centerconf[jj].getDouble("lat");
       centers.push_back(eckit::geometry::Point2(lon, lat));
       boost::shared_ptr<ObsSpaces_>
          lobs(new ObsSpaces_(obsdb, centers[jj], dist, max_nobs));
       localobs.push_back(lobs);
       Log::test() << "Local obs around: " << centers[jj] << std::endl;
       Log::test() << *localobs[jj] << std::endl;
       //  Setup obs bias<
       boost::shared_ptr<ObsAuxCtrls_> lobias(new ObsAuxCtrls_(obsdb, obsconf));
       localobias.push_back(lobias);
       //  Setup observer
       boost::shared_ptr<Observers<MODEL, State_>>
          lpobs(new Observers<MODEL, State_>(obsconf, *localobs[jj], *localobias[jj]));
       pobs.push_back(lpobs);
       post.enrollProcessor(pobs[jj]);
    }

//  Compute H(x)
    model.forecast(xx, moderr, winlen, post);
    Log::info() << "LocalHofX: Finished observation computation." << std::endl;
    Log::test() << "Final state: " << xx << std::endl;

//  Save H(x)
    for (std::size_t jj = 0; jj < centerconf.size(); ++jj) {
       const Observations_ & yobs = pobs[jj]->hofx();
       Log::test() << jj << " local H(x): " << yobs << std::endl;
       yobs.save("hofx");
    }
//  Read full H(x)
    Observations_ yobs(obsdb, "hofx");
    Log::test() << "H(x): " << yobs << std::endl;
    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::LocalHofX<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_LOCALHOFX_H_
