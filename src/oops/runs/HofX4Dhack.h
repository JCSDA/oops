/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_HOFX4DHACK_H_
#define OOPS_RUNS_HOFX4DHACK_H_

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/assimilation/CalcHofX.h"
#include "oops/base/Geometry.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/Variables.h"
#include "oops/generic/instantiateVariableChangeFactory.h"
#include "oops/interface/ChangeVariables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/GetValues.h"
#include "oops/interface/Locations.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

// ########################################################################
// ########################################################################
// ### This file is a temporary hack until all models can work with a   ###
// ### single variable transform in GetValues. Then it will be removed. ###
// ########################################################################
// ########################################################################

namespace oops {

template <typename MODEL, typename OBS> class HofX4Dhack : public Application {
  typedef CalcHofX<OBS>              CalcHofX_;
  typedef Geometry<MODEL>            Geometry_;
  typedef GeoVaLs<OBS>               GeoVaLs_;
  typedef GetValues<MODEL, OBS>      GetValues_;
  typedef Locations<OBS>             Locations_;
  typedef ObsAuxControls<OBS>        ObsAux_;
  typedef Observations<OBS>          Observations_;
  typedef ObsSpaces<OBS>             ObsSpaces_;
  typedef State<MODEL>               State_;
  typedef ChangeVariables<MODEL>     ChangeVariables_;

  typedef std::vector<std::unique_ptr<GeoVaLs_>>    GeoVaLsVec_;
  typedef std::vector<std::unique_ptr<Locations_>>  LocationsVec_;
  typedef std::vector<Variables> VariablesVec_;

 public:
// -----------------------------------------------------------------------------
  explicit HofX4Dhack(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateObsFilterFactory<OBS>();
    instantiateVariableChangeFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofX4Dhack() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    // Setup observation window
    const util::Duration winlen(fullConfig.getString("window length"));
    const util::DateTime winbgn(fullConfig.getString("window begin"));
    const util::DateTime winend(winbgn + winlen);
    Log::info() << "Observation window from " << winbgn << " to " << winend << std::endl;

    // Get information about states
    std::vector<eckit::LocalConfiguration> statesConf;
    if (fullConfig.has("states")) {
      statesConf = fullConfig.getSubConfigurations("states");
    } else {
      if (fullConfig.has("state")) {
        statesConf[0] = eckit::LocalConfiguration(fullConfig, "state");
      } else {
        eckit::BadParameter("Error in background state(s) configuration");
      }
    }
    size_t nsubwin = statesConf.size();

    size_t ntasks = this->getComm().size();
    size_t myrank = this->getComm().rank();
    ASSERT(ntasks % nsubwin == 0);
    size_t ntaskpslot = ntasks / nsubwin;
    size_t mysubwin = myrank / ntaskpslot;
    Log::debug() << "ntasks = " << ntasks << " nsubwin = " << nsubwin << std::endl;

    // Define local sub-window
    util::Duration subWinLength = winlen;
    if (nsubwin > 1) subWinLength = winlen / (nsubwin - 1);
    Log::debug() << "Task " << mysubwin << " subWinLength " << subWinLength << std::endl;
    util::DateTime subWinTime  = winbgn + mysubwin * subWinLength;
    util::DateTime subWinBegin = subWinTime - subWinLength/2;
    util::DateTime subWinEnd   = subWinTime + subWinLength/2;
    if (mysubwin == 0) subWinBegin = subWinTime;
    if (mysubwin == nsubwin - 1) subWinEnd = subWinTime;
    ASSERT(subWinBegin >= winbgn);
    ASSERT(subWinEnd <= winend);
    Log::debug() << "Task " << mysubwin << " Obs times " << subWinTime
                 << ", " << subWinBegin << ", " << subWinEnd << std::endl;
    // Create a communicator for same sub-window, to be used for communications in space
    std::string sgeom = "comm_geom_" + std::to_string(mysubwin);
    char const *geomName = sgeom.c_str();
    eckit::mpi::Comm * commSpace = &this->getComm().split(mysubwin, geomName);
    // Create a communicator for same local area, to be used for communications in time
    size_t myarea = commSpace->rank();
    std::string stime = "comm_time_" + std::to_string(myarea);
    char const *timeName = stime.c_str();
    eckit::mpi::Comm * commTime = &this->getComm().split(myarea, timeName);
    ASSERT(commTime->size() == nsubwin);

    // Setup geometry
    const eckit::LocalConfiguration geometryConfig(fullConfig, "geometry");
    const Geometry_ geometry(geometryConfig, *commSpace);

    // Setup states for H(x)
    Log::info() << "State configuration is:" << statesConf[mysubwin] << std::endl;
    State_ xx(geometry, statesConf[mysubwin]);
    Log::test() << "State: " << xx << std::endl;
    Log::debug() << "Task " << mysubwin << " State time " << xx.validTime() << std::endl;
    ASSERT(xx.validTime() == subWinTime);

    // Setup observations
    const eckit::LocalConfiguration obsConfig(fullConfig, "observations");
    ObsSpaces_ obspaces(obsConfig, *commSpace, subWinBegin, subWinEnd, *commTime);
    ObsAux_ obsaux(obspaces, obsConfig);
    CalcHofX_ hofx(obspaces, obsConfig);
    hofx.initialize(obsaux);

    // fill in GeoVaLs
    GeoVaLsVec_ geovals;
    const LocationsVec_ & locations = hofx.locations();
    const VariablesVec_ & vars = hofx.requiredVars();
    Log::debug() << "HofX4Dhack: Required hofx size = " << hofx.requiredVars().size() << std::endl;

    Variables geovars;
    Log::debug() << "HofX4Dhack: Required vars size = " << vars.size() << std::endl;
    for (size_t jj = 0; jj < vars.size(); ++jj) {
      Log::debug() << "HofX4Dhack: Required vars:" << vars[jj] << std::endl;
      geovars += vars[jj];
    }
    Log::debug() << "HofX4Dhack: Required variables:" << geovars << std::endl;
    eckit::LocalConfiguration chvarconf;  // empty for now
    ChangeVariables_ chvar(chvarconf, geometry, xx.variables(), geovars);

    State_ zz(geometry, geovars, xx.validTime());
    chvar.changeVar(xx, zz);

    std::vector<eckit::LocalConfiguration> getValuesConfig =
      util::vectoriseAndFilter(obsConfig, "get values");

     // loop over all observation types
    for (size_t jj = 0; jj < obspaces.size(); ++jj) {
      GetValues_ getvals(geometry, *locations[jj], getValuesConfig[jj]);
      geovals.emplace_back(new GeoVaLs_(*locations[jj], vars[jj],
                                        geometry.variableSizes(vars[jj])));
      getvals.fillGeoVaLs(zz, subWinBegin, subWinEnd, *geovals[jj]);
    }

    // Compute H(x) on filled in geovals and run the filters
    Observations_ yobs = hofx.compute(geovals);
    hofx.saveQcFlags("EffectiveQC");
    hofx.saveObsErrors("EffectiveError");

    // Save H(x)
    Log::test() << "H(x): " << std::endl << yobs << "End H(x)" << std::endl;
    Log::info() << "H(x): " << std::endl << yobs << "End H(x)" << std::endl;
    yobs.save("hofx");
    obspaces.save();

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::HofX4Dhack<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HOFX4DHACK_H_
