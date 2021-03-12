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
#include "oops/base/IncrementEnsemble4D.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/StateEnsemble4D.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/Increment.h"
#include "oops/mpi/mpi.h"
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
  typedef IncrementEnsemble4D<MODEL>       IncrementEnsemble4D_;
  typedef Increment<MODEL>                 Increment_;
  typedef LocalEnsembleSolver<MODEL, OBS>  LocalSolver_;
  typedef ObsSpaces<OBS>                   ObsSpaces_;
  typedef Observations<OBS>                Observations_;
  typedef State4D<MODEL>                   State4D_;
  typedef StateEnsemble4D<MODEL>           StateEnsemble4D_;

 public:
// -----------------------------------------------------------------------------

  explicit LocalEnsembleDA(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
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
    const eckit::LocalConfiguration obsConfig(fullConfig, "observations");
    ObsSpaces_ obsdb(obsConfig, this->getComm(), winbgn, winend);
    Observations_ yobs(obsdb, "ObsValue");

    // Get background configurations
    const eckit::LocalConfiguration bgConfig(fullConfig, "background");
    // Get driver configuration
    const eckit::LocalConfiguration driverConfig(fullConfig, "driver");

    // Read all ensemble members
    StateEnsemble4D_ ens_xx(geometry, bgConfig);
    const size_t nens = ens_xx.size();
    const Variables statevars = ens_xx.variables();

    // set up solver
    std::unique_ptr<LocalSolver_> solver =
         LocalEnsembleSolverFactory<MODEL, OBS>::create(obsdb, geometry, fullConfig, nens);

    // test prints for the prior ensemble
    bool do_test_prints = driverConfig.getBool("do test prints", true);
    if (do_test_prints) {
      for (size_t jj = 0; jj < nens; ++jj) {
        Log::test() << "Initial state for member " << jj+1 << ":" << ens_xx[jj] << std::endl;
      }
    }

    // compute H(x)
    bool readFromDisk = driverConfig.getBool("read HX from disk", false);
    Observations_ yb_mean = solver->computeHofX(ens_xx, 0, readFromDisk);
    Log::test() << "H(x) ensemble background mean: " << std::endl << yb_mean << std::endl;

    Departures_ ombg(yobs - yb_mean);
    ombg.save("ombg");
    Log::test() << "background y - H(x): " << std::endl << ombg << std::endl;

    // quit early if running in observer-only mode
    bool observerOnly = driverConfig.getBool("run as observer only", false);
    if (observerOnly) {return 0;}

    // calculate background mean
    State4D_ bkg_mean = ens_xx.mean();
    if (do_test_prints) {
      Log::test() << "Background mean :" << bkg_mean << std::endl;
    }

    // calculate background ensemble perturbations
    IncrementEnsemble4D_ bkg_pert(ens_xx, bkg_mean, statevars);

    // initialize empty analysis perturbations
    IncrementEnsemble4D_ ana_pert(geometry, statevars, ens_xx[0].validTimes(), bkg_pert.size());

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

    // save the posterior mean and ensemble first (since they are needed for the next cycle)
    // save the posterior mean
    State4D_ ana_mean = ens_xx.mean();   // calculate analysis mean
    if (do_test_prints) {
      Log::test() << "Analysis mean :" << ana_mean << std::endl;
    }
    bool save_xamean = driverConfig.getBool("save posterior mean", false);
    if (save_xamean) {
      eckit::LocalConfiguration outConfig(fullConfig, "output");
      outConfig.set("member", 0);
      ana_mean.write(outConfig);
    }

    // save the posterior ensemble
    bool save_ens = driverConfig.getBool("save posterior ensemble", true);
    if (save_ens) {
      size_t mymember;
      for (size_t jj=0; jj < nens; ++jj) {
        mymember = jj+1;
        eckit::LocalConfiguration outConfig(fullConfig, "output");
        outConfig.set("member", mymember);
        ens_xx[jj].write(outConfig);
      }
    }

    // below is the diagnostic output -----------------------------
    // save the background mean
    bool save_xbmean = driverConfig.getBool("save prior mean", false);
    if (save_xbmean) {
      eckit::LocalConfiguration outConfig(fullConfig, "output mean prior");
      outConfig.set("member", 0);
      bkg_mean.write(outConfig);
    }

    // save the analysis increment
    bool save_mean_increment = driverConfig.getBool("save posterior mean increment", false);
    if (save_mean_increment) {
      eckit::LocalConfiguration outConfig(fullConfig, "output increment");
      outConfig.set("member", 0);
      for (size_t itime = 0; itime < ana_mean.size(); ++itime) {
        Increment_ ana_increment(ana_pert[0][itime], false);
        ana_increment.diff(ana_mean[itime], bkg_mean[itime]);
        ana_increment.write(outConfig);
        if (do_test_prints) {
          Log::test() << "Analysis mean increment :" << ana_increment << std::endl;
        }
      }
    }

    // save the prior variance
    bool save_variance = driverConfig.getBool("save prior variance", false);
    if (save_variance) {
      eckit::LocalConfiguration outConfig(fullConfig, "output variance prior");
      outConfig.set("member", 0);
      std::string strOut("Forecast variance :");
      saveVariance(outConfig, bkg_pert, do_test_prints, strOut);
    }

    // save the posterior variance
    save_variance = driverConfig.getBool("save posterior variance", false);
    if (save_variance) {
      eckit::LocalConfiguration outConfig(fullConfig, "output variance posterior");
      outConfig.set("member", 0);
      std::string strOut("Analysis variance :");
      saveVariance(outConfig, ana_pert, do_test_prints, strOut);
    }

    // posterior observer
    // note: if H(X) is read from file, it might have used different time slots for observation
    // than LETKF background/analysis perturbations.
    // hence one might not expect that oman and omaf are comparable
    bool do_posterior_observer = driverConfig.getBool("do posterior observer", true);
    if (do_posterior_observer) {
      Observations_ ya_mean = solver->computeHofX(ens_xx, 1, false);
      Log::test() << "H(x) ensemble analysis mean: " << std::endl << ya_mean << std::endl;

      // calculate analysis obs departures
      Departures_ oman(yobs - ya_mean);
      oman.save("oman");
      Log::test() << "analysis y - H(x): " << std::endl << oman << std::endl;

      // display overall background/analysis RMS stats
      Log::test() << "ombg RMS: " << ombg.rms() << std::endl
                << "oman RMS: " << oman.rms() << std::endl;
    }

    return 0;
  }

// -----------------------------------------------------------------------------

 private:
  std::string appname() const {
    return "oops::LocalEnsembleDA<" + MODEL::name() + ", " + OBS::name() + ">";
  }

  void saveVariance(const eckit::LocalConfiguration & outConfig, const IncrementEnsemble4D_ & perts,
                    const bool do_test_prints, const std::string & strOut) const {
    // save and optionaly print varaince of an IncrementEnsemble4D_ object
    size_t nens = perts.size();
    const double nc = 1.0/(static_cast<double>(nens) - 1.0);
    for (size_t itime = 0; itime < perts[0].size(); ++itime) {
      Increment_ var(perts[0][itime], false);
      for (size_t iens = 0; iens < nens; ++iens) {
        Increment_ tmp(perts[iens][itime], true);
        tmp.schur_product_with(perts[iens][itime]);
        var += tmp;
      }
      var *= nc;
      var.write(outConfig);
      if (do_test_prints) {
        Log::test() << strOut << var << std::endl;
      }
    }
  }

// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_LOCALENSEMBLEDA_H_
