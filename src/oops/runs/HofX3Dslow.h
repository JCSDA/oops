/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_HOFX3DSLOW_H_
#define OOPS_RUNS_HOFX3DSLOW_H_

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CalcHofX.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
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

namespace oops {

template <typename MODEL, typename OBS> class HofX3Dslow : public Application {
  typedef CalcHofX<OBS>              CalcHofX_;
  typedef Geometry<MODEL>            Geometry_;
  typedef GeoVaLs<OBS>               GeoVaLs_;
  typedef GetValues<MODEL, OBS>      GetValues_;
  typedef Locations<OBS>             Locations_;
  typedef ObsAuxControls<OBS>        ObsAux_;
  typedef Observations<OBS>          Observations_;
  typedef ObsSpaces<OBS>             ObsSpaces_;
  typedef State<MODEL>               State_;

  typedef std::vector<std::unique_ptr<GeoVaLs_>>    GeoVaLsVec_;
  typedef std::vector<std::unique_ptr<Locations_>>  LocationsVec_;
  typedef std::vector<Variables> VariablesVec_;


 public:
// -----------------------------------------------------------------------------
  explicit HofX3Dslow(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateObsFilterFactory<OBS>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofX3Dslow() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    // Setup observation window
    const util::Duration winlen(fullConfig.getString("window length"));
    const util::DateTime winbgn(fullConfig.getString("window begin"));
    const util::DateTime winend(winbgn + winlen);
    Log::info() << "Observation window from " << winbgn << " to " << winend << std::endl;

    // Setup geometry
    const eckit::LocalConfiguration geometryConfig(fullConfig, "geometry");
    const Geometry_ geometry(geometryConfig, this->getComm());

    // Setup states for H(x)
    const eckit::LocalConfiguration stateConfig(fullConfig, "state");
    Log::info() << "State configuration is:" << stateConfig << std::endl;
    State_ xx(geometry, stateConfig);
    Log::test() << "State: " << xx << std::endl;

    // Setup observations
    const eckit::LocalConfiguration obsConfig(fullConfig, "observations");
    ObsSpaces_ obspaces(obsConfig, this->getComm(), winbgn, winend);
    ObsAux_ obsaux(obspaces, obsConfig);
    CalcHofX_ hofx(obspaces, obsConfig);
    hofx.initialize(obsaux);

    // fill in GeoVaLs
    GeoVaLsVec_ geovals;
    const LocationsVec_ & locations = hofx.locations();
    const VariablesVec_ & vars = hofx.requiredVars();

    std::vector<eckit::LocalConfiguration> getValuesConfig =
      util::vectoriseAndFilter(obsConfig, "get values");

     // loop over all observation types
    for (size_t jj = 0; jj < obspaces.size(); ++jj) {
      GetValues_ getvals(geometry, *locations[jj], getValuesConfig[jj]);
      geovals.emplace_back(new GeoVaLs_(*locations[jj], vars[jj]));
      getvals.fillGeoVaLs(xx, winbgn, winend, *geovals[jj]);
    }

    // Compute H(x) on filled in geovals and run the filters
    Observations_ yobs = hofx.compute(geovals);
    hofx.saveQcFlags("EffectiveQC");
    hofx.saveObsErrors("EffectiveError");

    // Save H(x)
    Log::test() << "H(x): " << std::endl << yobs << "End H(x)" << std::endl;
    yobs.save("hofx");

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::HofX3Dslow<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HOFX3DSLOW_H_
