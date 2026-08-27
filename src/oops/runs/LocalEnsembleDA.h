/*
 * (C) Copyright 2019-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_LOCALENSEMBLEDA_H_
#define OOPS_RUNS_LOCALENSEMBLEDA_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/config/YAMLConfiguration.h"
#include "oops/assimilation/instantiateLocalEnsembleSolverFactory.h"
#include "oops/assimilation/LocalEnsembleSolver.h"
#include "oops/base/Departures.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/StateSet.h"
#include "oops/base/StateSetSaver.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/runs/Forecast.h"
#include "oops/util/ConfigHelpers.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/printRunStats.h"


namespace oops {

// -----------------------------------------------------------------------------
/// \brief Application for local ensemble data assimilation
template <typename MODEL, typename OBS> class LocalEnsembleDA : public Application {
  typedef Departures<OBS>                  Departures_;
  typedef Geometry<MODEL>                  Geometry_;
  typedef GeometryIterator<MODEL>          GeometryIterator_;
  typedef Increment<MODEL>                 Increment_;
  typedef IncrementSet<MODEL>              IncrementSet_;
  typedef LocalEnsembleSolver<MODEL, OBS>  LocalSolver_;
  typedef ObsSpaces<OBS>                   ObsSpaces_;
  typedef Observations<OBS>                Observations_;
  typedef Model<MODEL>                     Model_;
  typedef ModelAuxControl<MODEL>           ModelAux_;
  typedef StateSet<MODEL>                  StateSet_;
  typedef State<MODEL>                     State_;

 public:
// -----------------------------------------------------------------------------

  explicit LocalEnsembleDA(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateLocalEnsembleSolverFactory<MODEL, OBS>();
  }

// -----------------------------------------------------------------------------

  virtual ~LocalEnsembleDA() = default;

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const override {
    std::unique_ptr<Geometry_> geometry;

    // Instantiate ens_xx depending on whether we are running inline or not
    auto ens_xx = [&] {
      if (fullConfig.getBool("Run Inline", false)) {
        auto object = localizeEnsembleFC(fullConfig, geometry);
        return object;
      } else {
        geometry = std::make_unique<Geometry_>(eckit::LocalConfiguration(fullConfig, "geometry"),
                                               this->getComm());
        auto object = StateSet_(*geometry, eckit::LocalConfiguration(fullConfig,
                                                                            "background"));
        return object;
      }
    }();

    //  Setup observation window
    const util::TimeWindow timeWindow(fullConfig.getSubConfiguration("time window"));
    Log::info() << "Observation window: " << timeWindow << std::endl;

    // Get observations configuration
    const eckit::LocalConfiguration observationsConfig(fullConfig, "observations");
    eckit::LocalConfiguration obsConfig = observationsConfig.getSubConfiguration("observers");

    // if any of the obs. spaces uses Halo distribution it will need to know the geometry
    // of the local grid on this PE
    if (fullConfig.getBool("driver.update obs config with geometry info", true))
        updateConfigWithPatchGeometry(*geometry, obsConfig);

    // Re-assemble the "observations" section around the (possibly modified) observers
    // list, so that ObsSpaces also sees the settings applying to every obs space.
    eckit::LocalConfiguration obsSpacesConfig(observationsConfig);
    obsSpacesConfig.set("observers", obsConfig.getSubConfigurations());

    // Setup observations
    const eckit::mpi::Comm & time = oops::mpi::myself();
    ObsSpaces_ obsdb(obsSpacesConfig, this->getComm(), timeWindow, time);

    // Read all ensemble members and compute the ensemble mean
    const size_t nens = ens_xx.ens_size();
    const Variables statevars = ens_xx.variables();
    Variables incvars;
    if (fullConfig.has("increment variables")) {
      incvars += Variables(fullConfig, "increment variables");
    } else {
      incvars += statevars;
    }
    StateSet_ bkg_mean = ens_xx.ens_mean();
    // if control member is present use that instead of the ensemble mean
    if (fullConfig.getBool("driver.use control member", false)) {
      StateSet_ controlMember(*geometry, eckit::LocalConfiguration(fullConfig, "control member"));
      bkg_mean = controlMember;
    }

    util::printRunStats("LocalEnsembleDA before solver ctor");

    // set up solver
    std::unique_ptr<LocalSolver_> solver =
         LocalEnsembleSolverFactory<MODEL, OBS>::create(obsdb, *geometry, fullConfig,
                                                        nens, bkg_mean, incvars);

    // test prints for the prior ensemble
    bool do_test_prints = fullConfig.getBool("driver.do test prints", true);
    if (do_test_prints) {
      for (size_t jj = 0; jj < nens; ++jj) {
        Log::test() << "Initial state for member " << jj+1 << ":";
        for (size_t jt = 0; jt < ens_xx.time_size(); ++jt) {
          Log::test() << ens_xx(jt, jj) << std::endl;
        }
      }
    }

    util::printRunStats("LocalEnsembleDA before computeHofX");

    // compute H(x)
    Observations_ yobs(obsdb, "ObsValue");
    Observations_ yb_mean = solver->computeHofX(ens_xx, 0,
                              fullConfig.getBool("driver.read HX from disk", false));
    if (do_test_prints) {
       Log::test() << "H(x) ensemble background mean: " << std::endl << yb_mean << std::endl;
    }

    Departures_ ombg(yobs - yb_mean);
    ombg.save("ombg");
    if (do_test_prints) {
       Log::test() << "background y - H(x): " << std::endl << ombg << std::endl;
    }

    // quit early if running in observer-only mode
    if (fullConfig.getBool("driver.run as observer only", false)) {
      obsdb.save();
      return 0;
    }

    // print background mean
    if (do_test_prints) {
      Log::test() << "Background mean :" << bkg_mean << std::endl;
    }

    std::vector<int> members(nens);
    std::iota(members.begin(), members.end(), 0);

    // calculate background ensemble perturbations
    IncrementSet_ bkg_pert(ens_xx, bkg_mean, incvars, members);

    // initialize empty analysis perturbations
    IncrementSet_ ana_pert(*geometry, incvars, ens_xx.times(), ens_xx.commTime(),
                           members);

    // run the solver at each gridpoint
    Log::info() << "Beginning core local solver..." << std::endl;
    util::printRunStats("LocalEnsembleDA before solver", true);

    solver->measurementUpdate(bkg_pert, ana_pert);

    // wait all tasks to finish their solution, so the timing for functions below reports
    // time which truly used (not from mpi_wait(), as all tasks need to sync before write).
    oops::mpi::world().barrier();

    Log::info() << "Local solver completed." << std::endl;
    util::printRunStats("LocalEnsembleDA after solver", true);
    // calculate final analysis states
    if (incvars == statevars) {
      // Initialize each ensemble member state from the background mean (per time slot)
      for (size_t jj = 0; jj < nens; ++jj) {
        for (size_t itime = 0; itime < ens_xx.time_size(); ++itime) {
          ens_xx(itime, jj) = bkg_mean[itime];
        }
      }
      // Apply all analysis perturbations across the ensemble/time set
      ens_xx += ana_pert;
    } else {
      // When increment variables differ, apply net analysis perturbations set-wise
      // by first forming analysis increments relative to background perturbations.
      IncrementSet_ net_ana(ana_pert);
      net_ana -= bkg_pert;
      ens_xx += net_ana;
    }

    // save the posterior mean, ensemble, and ensemble of increments first
    // (since they are needed for the next cycle)

    // save the posterior ensemble increments
    if (fullConfig.getBool("driver.save posterior ensemble increments", false)) {
      if (!fullConfig.has("output ensemble increments")) {
        throw eckit::BadValue(
          "`save posterior ensemble increment` is set to true, but `output ensemble increments` "
          "configuration not found.");
      }
      eckit::LocalConfiguration output(fullConfig, "output ensemble increments");
      for (size_t jj = 0; jj < nens; ++jj) {
        util::setMember(output, jj+1);
        for (size_t itime = 0; itime < ana_pert.time_size(); ++itime) {
          Increment_ ana_increment(ana_pert(itime, jj), true);
          ana_increment -= bkg_pert(itime, jj);
          ana_increment.write(output);
        }
      }
    }

    // save the posterior mean
    StateSet_ ana_mean = ens_xx.ens_mean();   // calculate analysis mean
    if (do_test_prints) {
      Log::test() << "Analysis mean :" << ana_mean << std::endl;
    }
    if (fullConfig.getBool("driver.save posterior mean", false)) {
      if (!fullConfig.has("output")) {
        throw eckit::BadValue("`save posterior mean` is set to true, but `output` "
                              "configuration not found.");
      }
      eckit::LocalConfiguration outConfig(fullConfig, "output");
      outConfig.set("member", 0);
      ana_mean.write(outConfig);
    }

    // save the posterior ensemble
    if (fullConfig.getBool("driver.save posterior ensemble", true)) {
      if (!fullConfig.has("output")) {
        throw eckit::BadValue("`save posterior ensemble` is set to true, but `output` "
                              "configuration not found.");
      }
      eckit::LocalConfiguration outConfig(fullConfig, "output");
      ens_xx.write(outConfig);
    }

    // below is the diagnostic output -----------------------------
    // save the background mean
    if (fullConfig.getBool("driver.save prior mean", false)) {
      if (!fullConfig.has("output mean prior")) {
        throw eckit::BadValue("`save prior mean` is set to true, but `output mean prior` "
                              "configuration not found.");
      }
      eckit::LocalConfiguration outConfig(fullConfig, "output mean prior");
      outConfig.set("member", 0);
      bkg_mean.write(outConfig);
    }

    // save the analysis mean increment
    if (fullConfig.getBool("driver.save posterior mean increment", false)) {
      if (!fullConfig.has("output increment")) {
        throw eckit::BadValue("`save posterior mean increment` is set to true, but "
                              "`output increment` configuration not found.");
      }
      eckit::LocalConfiguration output(fullConfig, "output increment");
      util::setMember(output, 0);
      for (size_t itime = 0; itime < ana_mean.time_size(); ++itime) {
        Increment_ ana_increment(ana_pert(itime, 0), false);
        ana_increment.diff(ana_mean[itime], bkg_mean[itime]);
        ana_increment.write(output);
        if (do_test_prints) {
          Log::test() << "Analysis mean increment :" << ana_increment << std::endl;
        }
      }
    }

    // save the prior variance
    if (fullConfig.getBool("driver.save prior variance", false)) {
      if (!fullConfig.has("output variance prior")) {
        throw eckit::BadValue("`save prior variance` is set to true, but `output variance prior` "
                              "configuration not found.");
      }
      eckit::LocalConfiguration output(fullConfig, "output variance prior");
      util::setMember(output, 0);
      std::string strOut("Forecast variance :");
      saveVariance(output, bkg_pert, do_test_prints, strOut);
    }

    // save the posterior variance
    if (fullConfig.getBool("driver.save posterior variance", false)) {
      if (!fullConfig.has("output variance posterior")) {
        throw eckit::BadValue("`save posterior variance` is set to true, but "
                              "`output variance posterior` configuration not found.");
      }
      eckit::LocalConfiguration output(fullConfig, "output variance posterior");
      util::setMember(output, 0);
      std::string strOut("Analysis variance :");
      saveVariance(output, ana_pert, do_test_prints, strOut);
    }

    // posterior observer
    // note: if H(X) is read from file, it might have used different time slots for observation
    // than LETKF background/analysis perturbations.
    // hence one might not expect that oman and omaf are comparable
    if (fullConfig.getBool("driver.do posterior observer", true)) {
      // need to create a posterior solver that stores ana_mean internally.
      // This is needed if linear observer is used, because it is linearized arround this mean
      std::unique_ptr<LocalSolver_> posteriorSolver =
         LocalEnsembleSolverFactory<MODEL, OBS>::create(obsdb, *geometry, fullConfig,
                                                        nens, ana_mean, incvars);
      Observations_ ya_mean = posteriorSolver->computeHofX(ens_xx, 1, false);
      Log::test() << "H(x) ensemble analysis mean: " << std::endl << ya_mean << std::endl;

      // calculate analysis obs departures
      Departures_ oman(yobs - ya_mean);
      oman.save("oman");
      Log::test() << "analysis y - H(x): " << std::endl << oman << std::endl;

      // display overall background/analysis RMS stats
      Log::test() << "ombg RMS: " << ombg.rms() << std::endl
                << "oman RMS: " << oman.rms() << std::endl;
    }

    // Save the obsspace only if an hofx was calculated
    // (either prior and/or posterior)
    if (!fullConfig.getBool("driver.read HX from disk", false) ||
        fullConfig.getBool("driver.do posterior observer", true)) {
      obsdb.save();
    }

    return 0;
  }

// -----------------------------------------------------------------------------

  StateSet_ localizeEnsembleFC(const eckit::Configuration & fullConfig,
          std::unique_ptr<Geometry_> & DAgeometry) const {
  // This function creates a DA geometry that has the same resolution as the forecast geometry, but
  // is decomposed into patches that are N times smaller than the forecast geometry, where N is the
  // number of ensemble members. Note that the DAgeometry layout must be evenly divisible by the
  // forecast layout. (e.g. DA layout = 4,4, FC layout = 2,2, N = 4)
  // Also note that DA layout nx * ny = FC layout nx * ny * N
  // Next, this will run a set of ensemble forecasts (or read in previously computed forecasts)
  // then "localize" the State variables and return a StateSet object with all of the state
  // variables on the local ensemble of the StateSet held by the StateSet variable returned

    // Get the MPI partition

    eckit::LocalConfiguration inlineParams = fullConfig.getSubConfiguration("inline parameters");
    const std::vector<std::string> files = inlineParams.getStringVector("Forecast configuration");
    const int batchSize = inlineParams.getInt("forecast batch size");
    const std::string pattern = inlineParams.getString("output file pattern");

    const int nmembers = files.size();
    const int ntasks = this->getComm().size();
    const int mytask = this->getComm().rank();  // global rank
    const int tasks_per_member = ntasks / nmembers;
    // divide by blocks of tasks_per_member
    int mymember = mytask / tasks_per_member + 1;

    eckit::LocalConfiguration subconfig = fullConfig.getSubConfiguration("geometry");
    // the layout here needs to be nmembers * the layout for the forecast geometry
    DAgeometry = std::unique_ptr<Geometry_>(new Geometry_(subconfig, this->getComm() ));

    Log::info() << "Running " << nmembers << " EnsembleGETKFApplication members handled by "
                << ntasks << " total MPI tasks and "
                << tasks_per_member << " MPI tasks per member." << std::endl;

    ASSERT(ntasks%nmembers == 0);

    //  Create the communicator for each ensemble member, named comm_member_{i}:
    std::string commNameStr = "comm_member_" + std::to_string(mymember);
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commMember = this->getComm().split(mymember, commName);
    const int subrank = commMember.rank();

    //  Create the communicator for each decomposed patch of geometry
    std::string patchNameStr = "patch_member_" + std::to_string(subrank);
    char const *patchName = patchNameStr.c_str();
    eckit::mpi::Comm & patchMember = this->getComm().split(subrank, patchName);

    Log::info() << "size of patchMember/ENS comm is " << patchMember.size() << std::endl;
    //  Each member uses a different configuration:
    eckit::PathName confPath = files[mymember-1];
    eckit::YAMLConfiguration memberConf(confPath);
    eckit::LocalConfiguration fcstparams = eckit::LocalConfiguration(memberConf);

    const Geometry_ FCgeometry(fcstparams.getSubConfiguration("geometry"), commMember);
    Log::info() << "done with geometry" << std::endl;

    //  Setup times
    Log::info() << "setting up times" << std::endl;
    eckit::LocalConfiguration model = fcstparams.getSubConfiguration("model");
    const util::Duration tstep(model.getString("tstep"));
    eckit::LocalConfiguration ic = fcstparams.getSubConfiguration("initial condition");
    const util::DateTime bgndate(ic.getString("datetime"));
    const util::Duration fclength(fcstparams.getString("forecast length"));
    const util::DateTime enddate(bgndate + fclength);
    std::vector<util::DateTime> times;
    const Variables vars(ic, "state variables");

    // Don't save the initial state
    oops::mpi::world().barrier();
    for (util::DateTime ii=(bgndate+tstep); ii <= enddate; ii=ii+tstep) {
       Log::info() << "pushing back time " << ii << std::endl;
       times.push_back(ii);
    }
    oops::mpi::world().barrier();
    std::vector<int> ens;  // vector of ensemble indices (0-based)
    for (int m = 1; m <= nmembers; ++m) { ens.push_back(m); }
    oops::mpi::world().barrier();

    std::unique_ptr<StateSet_> ens_SS;
    PostProcessor<State_> post;  // Create the post processor where StateSet will be stored
    StateSetSaver<MODEL> *saver_ =
        new StateSetSaver<MODEL>(memberConf, FCgeometry, times, oops::mpi::myself(),
                    ens, patchMember);
    post.enrollProcessor(saver_);
    //  Each member uses a different configuration:
    for (int m = 1; m <= nmembers; ++m) {
      if ( m == mymember ) {
        Log::info() << "running on mymember = " << mymember  << " " << mytask << std::endl;
        executeForecast(FCgeometry, memberConf, post);
        Log::info() << "Done with ens execute\n";
      }
      if ( batchSize > 0 ) {  // don't divide by zero
        if (m % batchSize == 0) oops::mpi::world().barrier();
      }
    }
    oops::mpi::world().barrier();
    ens_SS = std::move(saver_->getStateSet());

    // just finished the forecast on FCgeometry that has N times bigger patches than global DAgeom
    // Pull the values from the local FCgeometry and put them into DAgeom
    StateSet_ localVec = ens_SS->transpose(this->getComm(), *DAgeometry, mymember);
    return(localVec);
  }

// -----------------------------------------------------------------------------

 private:
  std::string appname() const override {
    return "oops::LocalEnsembleDA<" + MODEL::name() + ", " + OBS::name() + ">";
  }

  void calculate_patchCenter(const Geometry_ & geometry, std::vector<double> & patchCenter) const {
    eckit::geometry::Point3 gptmp3;
    const double deg2rad = 3.14159265/180.0;

    // compute patch center.
    // Convert from spherical lat-lon coordinates to Cartesian x,y,z coordinate frame
    // Calculate the mean xyz position,
    // and then convert this mean position back to spherical lat-lon,
    // and use this mean lat-lon as patch center.
    double alat = 0.0;
    double alon = 0.0;
    double xmean = 0.0;
    double ymean = 0.0;
    double zmean = 0.0;
    int n = 0;
    for (GeometryIterator_ i = geometry.begin(); i != geometry.end(); ++i) {
      gptmp3 = *i;
      alon = gptmp3[0]*deg2rad;
      alat = gptmp3[1]*deg2rad;
      xmean += cos(alat)*cos(alon);
      ymean += cos(alat)*sin(alon);
      zmean += sin(alat);
      ++n;
    }
    xmean = xmean/static_cast<double>(n);
    ymean = ymean/static_cast<double>(n);
    zmean = zmean/static_cast<double>(n);

    double rmean = sqrt(xmean*xmean + ymean*ymean);
    patchCenter[0] = atan2(ymean, xmean)/deg2rad;
    patchCenter[1] = atan2(zmean, rmean)/deg2rad;
  }

  void updateConfigWithPatchGeometry(const Geometry_ & geometry,
                                     eckit::LocalConfiguration & obsConfig) const {
    std::vector<double> patchCenter(2, 0.0);
    double patchRadius = 0.0;

    eckit::geometry::Point2 gptmp;
    const double radius_earth = 6.371e6;

    // Calculate region's patch center.
    calculate_patchCenter(geometry, patchCenter);

    // compute radius
    eckit::geometry::Point2 center(patchCenter[0], patchCenter[1]);
    for (GeometryIterator_ i = geometry.begin(); i != geometry.end(); ++i) {
      gptmp[0] = (*i)[0];
      gptmp[1] = (*i)[1];
      double dist = eckit::geometry::Sphere::distance(radius_earth, center, gptmp);
      patchRadius = fmax(patchRadius, dist);
    }

    // update observations configs with information on patch center and radius
    std::vector<eckit::LocalConfiguration> obsConfigs = obsConfig.getSubConfigurations();

    if (obsConfigs.size() > 0) {
      for (auto & conf : obsConfigs) {
        conf.set("obs space.distribution.center", patchCenter);
        conf.set("obs space.distribution.radius", patchRadius);
      }

      eckit::LocalConfiguration tmp;
      tmp.set("observers", obsConfigs);
      obsConfig = tmp.getSubConfiguration("observers");
    } else {
      obsConfig.set("obs space.distribution.center", patchCenter);
      obsConfig.set("obs space.distribution.radius", patchRadius);
    }
  }

  void saveVariance(const eckit::LocalConfiguration & params, const IncrementSet_ & perts,
                    const bool do_test_prints, const std::string & strOut) const {
      // save and optionally print variance of an IncrementSet_ object
      IncrementSet_ stddev = perts.ens_stddev();
      for (size_t itime = 0; itime < perts.time_size(); ++itime) {
          Increment_ var = stddev(itime, 0);
          var.schur_product_with(var);
          // write to disk and do test prints
          var.write(params);
          if (do_test_prints) {
            Log::test() << strOut << var << std::endl;
          }
      }
  }

// -----------------------------------------------------------------------------

  void executeForecast(const Geometry_ & geometry,
      const eckit::Configuration & fullConfig,
      PostProcessor<State_> & post) const {
//  Setup Model
    Log::info() << "Forecast:setting up model" << std::endl;
    const Model_ model(geometry, eckit::LocalConfiguration(fullConfig, "model"));

//  Setup initial state
    State_ xx(geometry, fullConfig.getSubConfiguration("initial condition"));

//  Setup augmented state
    const ModelAux_ moderr(geometry, fullConfig.getSubConfiguration("model aux control"));

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
};

}  // namespace oops
#endif  // OOPS_RUNS_LOCALENSEMBLEDA_H_
