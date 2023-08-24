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
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/instantiateLocalEnsembleSolverFactory.h"
#include "oops/assimilation/LocalEnsembleSolver.h"
#include "oops/base/Departures.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementEnsemble4D.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State4D.h"
#include "oops/base/StateEnsemble4D.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/printRunStats.h"


namespace oops {

// -----------------------------------------------------------------------------
/// \brief Options controlling output and observer for LocalEnsembleDA application.
class LocalEnsembleDADriverParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(LocalEnsembleDADriverParameters, Parameters)

 public:
  Parameter<bool> updateObsConfig{"update obs config with geometry info",
                  "controls whether observations config needs to be updated with information "
                  "about geometry distribution",
                  true, this};
  Parameter<bool> doTestPrints{"do test prints",
                  "controls whether additional output is printed to test stream",
                  true, this};
  Parameter<bool> readHofX{"read HX from disk",
                  "controls whether H(x) is read or computed", false, this};
  Parameter<bool> runObsOnly{"run as observer only",
                  "controls whether only observer is run, or observer and solver",
                  false, this};
  Parameter<bool> savePostMean{"save posterior mean",
                  "controls whether posterior (analysis) ensemble mean is saved",
                  false, this};
  Parameter<bool> savePostEns{"save posterior ensemble",
                  "controls whether posterior (analysis) ensemble is saved",
                  true, this};
  Parameter<bool> savePriorMean{"save prior mean",
                  "controls whether prior (background) ensemble mean is saved",
                  false, this};
  Parameter<bool> savePostMeanInc{"save posterior mean increment",
                  "controls whether posterior (analysis) mean ensemble increment is saved",
                  false, this};
  Parameter<bool> savePostEnsInc{"save posterior ensemble increments",
                  "controls whether posterior (analysis) increments are saved",
                  false, this};
  Parameter<bool> savePriorVar{"save prior variance",
                  "controls whether prior (background) ensemble variance is saved",
                  false, this};
  Parameter<bool> savePostVar{"save posterior variance",
                  "controls whether posterior (analysis) ensemble variance is saved",
                  false, this};
  Parameter<bool> doPostObs{"do posterior observer",
                  "controls whether H(x) is computed for the posterior (analysis) ensemble",
                  true, this};
  Parameter<bool> useControlMember{"use control member",
                  "use control member to center prior ensemble instead of the prior ensemble mean",
                  false, this};
};

// -----------------------------------------------------------------------------
/// \brief Top-level options taken by the LocalEnsembleDA application.
template <typename MODEL>
class LocalEnsembleDAParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(LocalEnsembleDAParameters, ApplicationParameters)

  typedef Geometry<MODEL>  Geometry_;
  typedef Increment<MODEL> Increment_;

 public:
  typedef typename Geometry_::Parameters_       GeometryParameters_;
  typedef typename Increment_::WriteParameters_ IncrementWriteParameters_;

  /// Only observations taken at times lying in the (`window begin`, `window begin` + `window
  /// length`] interval will be included in observation spaces.
  RequiredParameter<util::DateTime> windowBegin{"window begin", this};
  RequiredParameter<util::Duration> windowLength{"window length", this};

  /// A list whose elements determine treatment of observations from individual observation spaces.
  /// Note: current code changes this section; it isn't trivial to define this as Parameters for
  /// now.
  RequiredParameter<eckit::LocalConfiguration> observations{"observations", this};

  RequiredParameter<GeometryParameters_> geometry{"geometry",
          "geometry used for all of the ensemble members and increments", this};

  Parameter<LocalEnsembleDADriverParameters> driver{"driver",
          "options controlling output and observer runs", {}, this};

  RequiredParameter<eckit::LocalConfiguration> background{"background",
          "ensemble of backgrounds", this};

  OptionalParameter<Variables> incvars{"increment variables",
          "analysis increment variables", this};

  RequiredParameter<eckit::LocalConfiguration> localEnsDA{"local ensemble DA",
          "local ensemble DA solver and its options", this};

  /// Note: these Parameters have to be present if driver.useControlMember==true
  OptionalParameter<eckit::LocalConfiguration> controlMember{"control member",
          "control member that can be used insteead of the ensemble mean", this};

  /// Note: these Parameters have to be present if driver.savePostMean or driver.savePostEns
  /// are true.
  OptionalParameter<eckit::LocalConfiguration> output{"output",
          "parameters for posterior mean and ensemble output", this};

  /// Note: these Parameters have to be present if driver.savePriorMean is true.
  OptionalParameter<eckit::LocalConfiguration> outputPriorMean{"output mean prior",
         "parameters for prior mean output", this};

  /// Note: these Parameters have to be present if driver.savePostMeanInc is true.
  OptionalParameter<IncrementWriteParameters_> outputPostMeanInc{"output increment",
         "parameters for posterior mean increment output", this};

  /// Note: these Parameters have to be present if driver.savePostEnsInc is true.
  OptionalParameter<IncrementWriteParameters_> outputPostEnsInc{"output ensemble increments",
         "parameters for posterior ensemble increments output", this};

  /// Note: these Parameters have to be present if driver.savePriorVar is true.
  OptionalParameter<IncrementWriteParameters_> outputPriorVar{"output variance prior",
         "parameters for prior variance output", this};

  /// Note: these Parameters have to be present if driver.savePostVar is true.
  OptionalParameter<IncrementWriteParameters_> outputPostVar{"output variance posterior",
         "parameters for posterior variance output", this};
};

// -----------------------------------------------------------------------------
/// \brief Application for local ensemble data assimilation
template <typename MODEL, typename OBS> class LocalEnsembleDA : public Application {
  typedef Departures<OBS>                  Departures_;
  typedef Geometry<MODEL>                  Geometry_;
  typedef GeometryIterator<MODEL>          GeometryIterator_;
  typedef IncrementEnsemble4D<MODEL>       IncrementEnsemble4D_;
  typedef Increment<MODEL>                 Increment_;
  typedef Increment4D<MODEL>               Increment4D_;
  typedef LocalEnsembleSolver<MODEL, OBS>  LocalSolver_;
  typedef ObsSpaces<OBS>                   ObsSpaces_;
  typedef Observations<OBS>                Observations_;
  typedef State4D<MODEL>                   State4D_;
  typedef StateEnsemble4D<MODEL>           StateEnsemble4D_;
  typedef typename Increment<MODEL>::WriteParameters_ IncrementWriteParameters_;
  typedef LocalEnsembleDAParameters<MODEL> LocalEnsembleDAParameters_;

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

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    // Deserialize parameters
    LocalEnsembleDAParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    //  Setup observation window
    const util::DateTime winbgn = params.windowBegin;
    const util::Duration winlen = params.windowLength;
    const util::DateTime winend(winbgn + winlen);
    Log::info() << "Observation window from " << winbgn << " to " << winend << std::endl;

    // Setup geometry
    const Geometry_ geometry(params.geometry, this->getComm());

    // Get observations configuration
    const eckit::LocalConfiguration observationsConfig = params.observations;
    eckit::LocalConfiguration obsConfig = observationsConfig.getSubConfiguration("observers");

    // if any of the obs. spaces uses Halo distribution it will need to know the geometry
    // of the local grid on this PE
    if (params.driver.value().updateObsConfig) updateConfigWithPatchGeometry(geometry, obsConfig);

    // Setup observations
    const eckit::mpi::Comm & time = oops::mpi::myself();
    ObsSpaces_ obsdb(obsConfig, this->getComm(), winbgn, winend, time);
    Observations_ yobs(obsdb, "ObsValue");

    // Read all ensemble members and compute the ensemble mean
    StateEnsemble4D_ ens_xx(geometry, params.background);
    const size_t nens = ens_xx.size();
    const Variables statevars = ens_xx.variables();
    Variables incvars;
    if (params.incvars.value() == boost::none) {
      incvars += statevars;
    } else {
      incvars += *params.incvars.value();
    }
    State4D_ bkg_mean = ens_xx.mean();
    // if control member is present use that instead of the ensemble mean
    if (params.driver.value().useControlMember) {
      State4D_ controlMember(geometry, *params.controlMember.value());
      bkg_mean = controlMember;
    }

    util::printRunStats("LocalEnsembleDA before solver ctor");

    // set up solver
    std::unique_ptr<LocalSolver_> solver =
         LocalEnsembleSolverFactory<MODEL, OBS>::create(obsdb, geometry, fullConfig,
                                                        nens, bkg_mean, incvars);

    // test prints for the prior ensemble
    bool do_test_prints = params.driver.value().doTestPrints;
    if (do_test_prints) {
      for (size_t jj = 0; jj < nens; ++jj) {
        Log::test() << "Initial state for member " << jj+1 << ":" << ens_xx[jj] << std::endl;
      }
    }

    util::printRunStats("LocalEnsembleDA before computeHofX");

    // compute H(x)
    Observations_ yb_mean = solver->computeHofX(ens_xx, 0, params.driver.value().readHofX);
    if (do_test_prints) {
       Log::test() << "H(x) ensemble background mean: " << std::endl << yb_mean << std::endl;
    }

    Departures_ ombg(yobs - yb_mean);
    ombg.save("ombg");
    if (do_test_prints) {
       Log::test() << "background y - H(x): " << std::endl << ombg << std::endl;
    }

    // quit early if running in observer-only mode
    if (params.driver.value().runObsOnly.value()) {
      obsdb.save();
      return 0;
    }

    // print background mean
    if (do_test_prints) {
      Log::test() << "Background mean :" << bkg_mean << std::endl;
    }

    // calculate background ensemble perturbations
    IncrementEnsemble4D_ bkg_pert(ens_xx, bkg_mean, incvars);

    // initialize empty analysis perturbations
    IncrementEnsemble4D_ ana_pert(geometry, incvars, ens_xx[0].validTimes(), bkg_pert.size());

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
      for (size_t jj = 0; jj < nens; ++jj) {
        ens_xx[jj] = bkg_mean;
        ens_xx[jj] += ana_pert[jj];
      }
    } else {
      Increment4D_ ana_increment(geometry, incvars, ens_xx[0].validTimes());
      for (size_t jj = 0; jj < nens; ++jj) {
        ana_increment = ana_pert[jj];
        for (size_t itime = 0; itime < bkg_pert[jj].size(); ++itime) {
          ana_increment[itime] -= bkg_pert[jj][itime];
        }
        ens_xx[jj] += ana_increment;
      }
    }

    // save the posterior mean, ensemble, and ensemble of increments first
    // (since they are needed for the next cycle)

    // save the posterior ensemble increments
    if (params.driver.value().savePostEnsInc.value()) {
      if (params.outputPostEnsInc.value() == boost::none) {
        throw eckit::BadValue(
          "`save posterior ensemble increment` is set to true, but `output ensemble increments` "
          "configuration not found.");
      }
      IncrementWriteParameters_ output = *params.outputPostEnsInc.value();
      for (size_t jj = 0; jj < nens; ++jj) {
        output.setMember(jj+1);
        for (size_t itime = 0; itime < ana_pert[0].size(); ++itime) {
          Increment_ ana_increment(ana_pert[jj][itime], true);
          ana_increment -= bkg_pert[jj][itime];
          ana_increment.write(output);
        }
      }
    }

    // save the posterior mean
    State4D_ ana_mean = ens_xx.mean();   // calculate analysis mean
    if (do_test_prints) {
      Log::test() << "Analysis mean :" << ana_mean << std::endl;
    }
    if (params.driver.value().savePostMean.value()) {
      if (params.output.value() == boost::none) {
        throw eckit::BadValue("`save posterior mean` is set to true, but `output` "
                              "configuration not found.");
      }
      eckit::LocalConfiguration outConfig = *params.output.value();
      outConfig.set("member", 0);
      ana_mean.write(outConfig);
    }

    // save the posterior ensemble
    if (params.driver.value().savePostEns.value()) {
      if (params.output.value() == boost::none) {
        throw eckit::BadValue("`save posterior ensemble` is set to true, but `output` "
                              "configuration not found.");
      }
      eckit::LocalConfiguration outConfig = *params.output.value();
      for (size_t jj = 0; jj < nens; ++jj) {
        outConfig.set("member", jj+1);
        ens_xx[jj].write(outConfig);
      }
    }

    // below is the diagnostic output -----------------------------
    // save the background mean
    if (params.driver.value().savePriorMean.value()) {
      if (params.outputPriorMean.value() == boost::none) {
        throw eckit::BadValue("`save prior mean` is set to true, but `output mean prior` "
                              "configuration not found.");
      }
      eckit::LocalConfiguration outConfig = *params.outputPriorMean.value();
      outConfig.set("member", 0);
      bkg_mean.write(outConfig);
    }

    // save the analysis mean increment
    if (params.driver.value().savePostMeanInc.value()) {
      if (params.outputPostMeanInc.value() == boost::none) {
        throw eckit::BadValue("`save posterior mean increment` is set to true, but "
                              "`output increment` configuration not found.");
      }
      IncrementWriteParameters_ output = *params.outputPostMeanInc.value();
      output.setMember(0);
      for (size_t itime = 0; itime < ana_mean.size(); ++itime) {
        Increment_ ana_increment(ana_pert[0][itime], false);
        ana_increment.diff(ana_mean[itime], bkg_mean[itime]);
        ana_increment.write(output);
        if (do_test_prints) {
          Log::test() << "Analysis mean increment :" << ana_increment << std::endl;
        }
      }
    }

    // save the prior variance
    if (params.driver.value().savePriorVar.value()) {
      if (params.outputPriorVar.value() == boost::none) {
        throw eckit::BadValue("`save prior variance` is set to true, but `output variance prior` "
                              "configuration not found.");
      }
      IncrementWriteParameters_ output = *params.outputPriorVar.value();
      output.setMember(0);
      std::string strOut("Forecast variance :");
      saveVariance(output, bkg_pert, do_test_prints, strOut);
    }

    // save the posterior variance
    if (params.driver.value().savePostVar.value()) {
      if (params.outputPostVar.value() == boost::none) {
        throw eckit::BadValue("`save posterior variance` is set to true, but "
                              "`output variance posterior` configuration not found.");
      }
      IncrementWriteParameters_ output = *params.outputPostVar.value();
      output.setMember(0);
      std::string strOut("Analysis variance :");
      saveVariance(output, ana_pert, do_test_prints, strOut);
    }

    // posterior observer
    // note: if H(X) is read from file, it might have used different time slots for observation
    // than LETKF background/analysis perturbations.
    // hence one might not expect that oman and omaf are comparable
    if (params.driver.value().doPostObs.value()) {
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

    // Save the obsspace only if an hofx was calculated
    // (either prior and/or posterior)
    if ( !params.driver.value().readHofX.value() ||
         params.driver.value().doPostObs.value()) {
      obsdb.save();
    }

    return 0;
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

  void saveVariance(const IncrementWriteParameters_ & params, const IncrementEnsemble4D_ & perts,
                    const bool do_test_prints, const std::string & strOut) const {
    // save and optionaly print varaince of an IncrementEnsemble4D_ object
    size_t nens = perts.size();
    const double ncVar = 1.0/(static_cast<double>(nens) - 1.0);
    const double ncMean = 1.0/(static_cast<double>(nens));
    for (size_t itime = 0; itime < perts[0].size(); ++itime) {
      // compute the mean
      Increment_ mean(perts[0][itime], false);
      for (size_t iens = 0; iens < nens; ++iens) {
         mean += perts[iens][itime];
      }
      mean *= ncMean;
      // remove the mean from the ensemble and accumulate sum of squares
      Increment_ var(perts[0][itime], false);
      for (size_t iens = 0; iens < nens; ++iens) {
        Increment_ tmp(perts[iens][itime], true);
        tmp -= mean;
        tmp.schur_product_with(tmp);
        var += tmp;
      }
      var *= ncVar;
      // write to disk and do test prints
      var.write(params);
      if (do_test_prints) {
        Log::test() << strOut << var << std::endl;
      }
    }
  }

// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    LocalEnsembleDAParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    LocalEnsembleDAParameters_ params;
    params.validate(fullConfig);
  }

// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_LOCALENSEMBLEDA_H_
