/*
 * (C) Copyright 2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_LETKF_H_
#define OOPS_RUNS_LETKF_H_

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CalcHofX.h"
#include "oops/assimilation/instantiateLETKFSolverFactory.h"
#include "oops/assimilation/LETKFSolverBase.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/StateEnsemble.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"


namespace oops {

/// Local Ensemble Tranform Kalman Filter (LETKF)
/*!
 * An (in progress) implementation of the LETKF from Hunt et al. 2007
 *
 * Hunt, B. R., Kostelich, E. J., & Szunyogh, I. (2007). Efficient data
 * assimilation for spatiotemporal chaos: A local ensemble transform Kalman
 * filter. Physica D: Nonlinear Phenomena, 230(1-2), 112-126.
 */
template <typename MODEL, typename OBS> class LETKF : public Application {
  typedef Departures<OBS>              Departures_;
  typedef DeparturesEnsemble<OBS>      DeparturesEnsemble_;
  typedef Geometry<MODEL>              Geometry_;
  typedef GeometryIterator<MODEL>      GeometryIterator_;
  typedef Increment<MODEL>             Increment_;
  typedef IncrementEnsemble<MODEL>     IncrementEnsemble_;
  typedef LETKFSolverBase<MODEL, OBS>  LETKFSolver_;
  typedef ObsEnsemble<OBS>             ObsEnsemble_;
  typedef ObsErrors<OBS>               ObsErrors_;
  typedef ObsSpaces<OBS>               ObsSpaces_;
  typedef Observations<OBS>            Observations_;
  typedef State<MODEL>                 State_;
  typedef State4D<MODEL>               State4D_;
  typedef StateEnsemble<MODEL>         StateEnsemble_;

 public:
// -----------------------------------------------------------------------------

  explicit LETKF(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateLETKFSolverFactory<MODEL, OBS>();
    instantiateObsErrorFactory<OBS>();
    instantiateObsFilterFactory<OBS>();
  }

// -----------------------------------------------------------------------------

  virtual ~LETKF() {}

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const {
    // Setup observation window
    const eckit::LocalConfiguration windowConfig(fullConfig, "Assimilation Window");
    const util::Duration winlen(windowConfig.getString("window_length"));
    const util::DateTime winbgn(windowConfig.getString("window_begin"));
    const util::DateTime winend(winbgn + winlen);
    const util::DateTime winhalf = winbgn + winlen/2;
    Log::debug() << "Observation window is: " << windowConfig << std::endl;

    // Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "Geometry");
    const Geometry_ resol(resolConfig, this->getComm());

    // Setup observations
    const eckit::LocalConfiguration obsConfig(fullConfig, "Observations");
    Log::debug() << "Observation configuration is: " << obsConfig << std::endl;
    ObsSpaces_ obsdb(obsConfig, this->getComm(), winbgn, winend);
    Observations_ yobs(obsdb, "ObsValue");

    // Get background configurations
    const eckit::LocalConfiguration bgConfig(fullConfig, "Background");
    const Variables statevars(bgConfig);

    // Read all ensemble members
    StateEnsemble_ ens_xx(resol, statevars, bgConfig);
    const size_t nens = ens_xx.size();
    ObsEnsemble_ obsens(obsdb, nens);

    // Initialize observer
    CalcHofX<MODEL, OBS> hofx(obsdb, resol, fullConfig);
    for (size_t jj = 0; jj < nens; ++jj) {
      // TODO(Travis) change the way input file name is specified, make
      //  more similar to how the output ens config is done
      Log::test() << "Initial state for member " << jj+1 << ":" << ens_xx[jj] << std::endl;

      // compute and save H(x)
      obsens[jj] = hofx.compute(ens_xx[jj]);
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << obsens[jj] << std::endl;
      obsens[jj].save("hofx0_"+std::to_string(jj+1));
    }
    // TODO(someone) still need to use QC flags (mask obsens)
    // QC flags and Obs errors are set to that of the last
    // ensemble member (those obs errors will be used in the assimilation)
    for (size_t jobs = 0; jobs < obsdb.size(); ++jobs) {
      hofx.qcFlags(jobs).save("EffectiveQC");
      hofx.obsErrors(jobs).save("EffectiveError");
    }

    // calculate background mean
    State4D_ bkg_mean = ens_xx.mean();
    Log::test() << "Background mean :" << bkg_mean << std::endl;

    // calculate background ensemble perturbations
    IncrementEnsemble_ bkg_pert(ens_xx, bkg_mean, statevars);

    // TODO(Travis) optionally save the background mean / standard deviation

    // calculate H(x) ensemble mean
    Observations_ yb_mean(obsens.mean());
    Log::test() << "H(x) ensemble background mean: " << std::endl << yb_mean << std::endl;

    // calculate H(x) ensemble perturbations
    DeparturesEnsemble_ ens_Yb(obsdb, nens);
    for (size_t iens = 0; iens < nens; ++iens) {
      ens_Yb[iens] = obsens[iens] - yb_mean;
    }

    // calculate obs departures
    Departures_ ombg(yobs - yb_mean);
    ombg.save("ombg");
    Log::test() << "background y - H(x): " << std::endl << ombg << std::endl;

    // initialize empty analysis perturbations
    IncrementEnsemble_ ana_pert(resol, statevars, ens_xx[0].validTimes(), bkg_pert.size());

    // set up solver
    std::unique_ptr<LETKFSolver_> localSolver =
         LETKFSolverFactory<MODEL, OBS>::create(fullConfig, nens);

    // run the LETKF solver at each gridpoint
    Log::info() << "Beginning core LETKF solver..." << std::endl;
    for (GeometryIterator_ i = resol.begin(); i != resol.end(); ++i) {
      // create the local subset of observations
      ObsSpaces_ local_obs(obsdb, *i, obsConfig);
      Departures_ local_ombg(local_obs, ombg);
      DeparturesEnsemble_ local_ens_Yb(local_obs, ens_Yb);
      // create local obs errors
      ObsErrors_ local_rmat(obsConfig, local_obs);

      localSolver->measurementUpdate(local_ombg, local_ens_Yb, local_rmat,
                                    bkg_pert, i, ana_pert);
    }
    Log::info() << "LETKF solver completed." << std::endl;

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

    // calculate oman
    for (size_t jj = 0; jj < nens; ++jj) {
      // compute and save H(x)
      obsens[jj] = hofx.compute(ens_xx[jj]);
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << obsens[jj] << std::endl;
      obsens[jj].save("hofx1_"+std::to_string(jj+1));
    }

    // calcualate H(x) ensemble analysis mean
    Observations_ ya_mean(obsens.mean());
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
    return "oops::LETKF<" + MODEL::name() + ", " + OBS::name() + ">";
  }

// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_LETKF_H_
