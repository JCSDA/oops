/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_HOFXNOMODEL_H_
#define OOPS_RUNS_HOFXNOMODEL_H_

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CalcHofX.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/GetValues.h"
#include "oops/interface/Locations.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL, typename OBS> class HofXNoModel : public Application {
  typedef CalcHofX<OBS>              CalcHofX_;
  typedef Geometry<MODEL>            Geometry_;
  typedef GeoVaLs<OBS>               GeoVaLs_;
  typedef GetValues<MODEL, OBS>      GetValues_;
  typedef Locations<OBS>             Locations_;
  typedef ObsAuxControls<OBS>        ObsAux_;
  typedef Observations<OBS>          Observations_;
  typedef ObsSpaces<OBS>             ObsSpaces_;
  typedef State4D<MODEL>             State4D_;

  typedef std::vector<std::unique_ptr<GeoVaLs_>>    GeoVaLsVec_;
  typedef std::vector<std::unique_ptr<Locations_>>  LocationsVec_;
  typedef std::vector<Variables> VariablesVec_;


 public:
// -----------------------------------------------------------------------------
  explicit HofXNoModel(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateObsFilterFactory<OBS>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofXNoModel() {}
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
    const eckit::LocalConfiguration stateConfig(fullConfig, "forecasts");
    Log::info() << "States configuration is:" << stateConfig << std::endl;
    State4D_ xx(geometry, stateConfig);
    Log::test() << "Initial state: " << xx[0] << std::endl;

    // Setup observations
    const eckit::LocalConfiguration obsConfig(fullConfig, "observations");
    ObsSpaces_ obspaces(obsConfig, this->getComm(), winbgn, winend);
    ObsAux_ obsaux(obspaces, obsConfig);
    CalcHofX_ hofx(obspaces, obsConfig);
    hofx.initialize(obsaux);

    // Setup and check time steps between the states in 4D state
    const size_t nstates = xx.size();
    util::Duration tstep = winlen;      // for a single state
    // if using several states, compute the timestep and check that it's the same for all states
    if (nstates > 1) {
      tstep = xx[1].validTime() - xx[0].validTime();
      for (size_t ii = 1; ii < nstates; ++ii) {
        ASSERT(tstep == (xx[ii].validTime() - xx[ii-1].validTime()));
      }
    }

    // fill in GeoVaLs
    GeoVaLsVec_ geovals;
    const LocationsVec_ & locations = hofx.locations();
    const VariablesVec_ & vars = hofx.requiredVars();

    std::vector<eckit::LocalConfiguration> getValuesConfig =
      util::vectoriseAndFilter(obsConfig, "get values");

     // loop over all observation types
    for (size_t jj = 0; jj < obspaces.size(); ++jj) {
      GetValues_ getvals(geometry, *locations[jj], getValuesConfig[jj]);
      // add GeoVaLs for this obs type
      geovals.emplace_back(new GeoVaLs_(*locations[jj], vars[jj]));
      // fill in by looping through 4D state
      for (size_t ii = 0; ii < nstates; ++ii) {
        util::DateTime t1 = std::max(xx[ii].validTime()-tstep/2, winbgn);
        util::DateTime t2 = std::min(xx[ii].validTime()+tstep/2, winend);
        getvals.fillGeoVaLs(xx[ii], t1, t2, *geovals[jj]);
      }
    }

    // Compute H(x) on filled in geovals and run the filters
    Observations_ yobs = hofx.compute(geovals);
    hofx.saveQcFlags("EffectiveQC");
    hofx.saveObsErrors("EffectiveError");
    Log::test() << "Final state: " << xx[xx.size()-1] << std::endl;

    // Save H(x)
    Log::test() << "H(x): " << std::endl << yobs << "End H(x)" << std::endl;
    yobs.save("hofx");

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::HofXNoModel<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HOFXNOMODEL_H_
