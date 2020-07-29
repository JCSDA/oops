/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_LOCALENSEMBLEDA_H_
#define OOPS_RUNS_LOCALENSEMBLEDA_H_

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/instantiateLocalEnsembleSolverFactory.h"
#include "oops/assimilation/LocalEnsembleSolver.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/Departures.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/StateEnsemble.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"


namespace oops {

/// \brief Application for local ensemble data assimilation
template <typename MODEL, typename OBS> class LocalEnsembleDA : public Application {
  typedef Departures<OBS>                  Departures_;
  typedef Geometry<MODEL>                  Geometry_;
  typedef GeometryIterator<MODEL>          GeometryIterator_;
  typedef IncrementEnsemble<MODEL>         IncrementEnsemble_;
  typedef LocalEnsembleSolver<MODEL, OBS>  LocalSolver_;
  typedef ObsSpaces<OBS>                   ObsSpaces_;
  typedef Observations<OBS>                Observations_;
  typedef State4D<MODEL>                   State4D_;
  typedef StateEnsemble<MODEL>             StateEnsemble_;

 public:
// -----------------------------------------------------------------------------

  explicit LocalEnsembleDA(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateLocalEnsembleSolverFactory<MODEL, OBS>();
    instantiateObsErrorFactory<OBS>();
    instantiateObsFilterFactory<OBS>();
  }

// -----------------------------------------------------------------------------

  virtual ~LocalEnsembleDA() = default;

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const {
    //  Setup observation window
    const util::DateTime winbgn(fullConfig.getString("window begin"));
    const util::Duration winlen(fullConfig.getString("window length"));
    const util::DateTime winend(winbgn + winlen);
    const util::DateTime winhalf = winbgn + winlen/2;
    Log::info() << "Observation window from " << winbgn << " to " << winend << std::endl;

    // Setup geometry
    const eckit::LocalConfiguration geometryConfig(fullConfig, "geometry");
    const Geometry_ geometry(geometryConfig, this->getComm());

    // Setup observations
    ObsSpaces_ obsdb(fullConfig, this->getComm(), winbgn, winend);
    Observations_ yobs(obsdb, "ObsValue");

    // Get background configurations
    const eckit::LocalConfiguration bgConfig(fullConfig, "background");

    // Read all ensemble members
    StateEnsemble_ ens_xx(geometry, bgConfig);
    const size_t nens = ens_xx.size();

    const Variables statevars = ens_xx.variables();

    // set up solver
    std::unique_ptr<LocalSolver_> solver =
         LocalEnsembleSolverFactory<MODEL, OBS>::create(obsdb, geometry, fullConfig, nens);

    for (size_t jj = 0; jj < nens; ++jj) {
      // TODO(Travis) change the way input file name is specified, make
      //  more similar to how the output ens config is done
      Log::test() << "Initial state for member " << jj+1 << ":" << ens_xx[jj] << std::endl;
    }

    // calculate background mean
    State4D_ bkg_mean = ens_xx.mean();
    Log::test() << "Background mean :" << bkg_mean << std::endl;

    // calculate background ensemble perturbations
    IncrementEnsemble_ bkg_pert(ens_xx, bkg_mean, statevars);

    // TODO(Travis) optionally save the background mean / standard deviation

    // compute H(x)
    Observations_ yb_mean = solver->computeHofX(ens_xx, 0);
    Log::test() << "H(x) ensemble background mean: " << std::endl << yb_mean << std::endl;

    Departures_ ombg(yobs - yb_mean);
    ombg.save("ombg");
    Log::test() << "background y - H(x): " << std::endl << ombg << std::endl;

    // initialize empty analysis perturbations
    IncrementEnsemble_ ana_pert(geometry, statevars, ens_xx[0].validTimes(), bkg_pert.size());

    // run the solver at each gridpoint
    Log::info() << "Beginning core local solver..." << std::endl;
    for (GeometryIterator_ i = geometry.begin(); i != geometry.end(); ++i) {
      solver->measurementUpdate(bkg_pert, i, ana_pert);
    }
    Log::info() << "Local solver completed." << std::endl;

    // calculate final analysis states
    for (size_t jj = 0; jj < nens; ++jj) {
      ens_xx[jj] = bkg_mean;
      ens_xx[jj] += ana_pert[jj];
    }

    // TODO(Travis) optionally save analysis standard deviation

    // save the analysis mean
    State4D_ ana_mean = ens_xx.mean();   // calculate analysis mean
    Log::info() << "Analysis mean :" << ana_mean << std::endl;
    eckit::LocalConfiguration outConfig(fullConfig, "output");
    outConfig.set("member", 0);
    ana_mean.write(outConfig);

    // save the analysis ensemble
    size_t mymember;
    for (size_t jj=0; jj < nens; ++jj) {
      mymember = jj+1;
      eckit::LocalConfiguration outConfig(fullConfig, "output");
      outConfig.set("member", mymember);
      ens_xx[jj].write(outConfig);
    }

    Observations_ ya_mean = solver->computeHofX(ens_xx, 1);
    Log::test() << "H(x) ensemble analysis mean: " << std::endl << ya_mean << std::endl;

    // calculate analysis obs departures
    Departures_ oman(yobs - ya_mean);
    oman.save("oman");
    Log::test() << "analysis y - H(x): " << std::endl << oman << std::endl;

    // display overall background/analysis RMS stats
    Log::test() << "ombg RMS: " << ombg.rms() << std::endl
                << "oman RMS: " << oman.rms() << std::endl;

    return 0;
  }

// -----------------------------------------------------------------------------

 private:
  std::string appname() const {
    return "oops::LocalEnsembleDA<" + MODEL::name() + ", " + OBS::name() + ">";
  }

// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_LOCALENSEMBLEDA_H_
