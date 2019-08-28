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

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Departures.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
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
template <typename MODEL> class LETKF : public Application {
  typedef ModelAuxControl<MODEL>  ModelAux_;
  typedef Departures<MODEL>        Departures_;
  typedef Geometry<MODEL>          Geometry_;
  typedef Increment<MODEL>         Increment_;
  typedef Model<MODEL>             Model_;
  typedef ObsAuxControls<MODEL>    ObsAuxCtrls_;
  typedef ObsEnsemble<MODEL>       ObsEnsemble_;
  typedef ObsErrors<MODEL>         ObsErrors_;
  typedef ObsSpaces<MODEL>         ObsSpaces_;
  typedef Observations<MODEL>      Observations_;
  typedef State<MODEL>             State_;
  template <typename DATA> using ObsData_ = ObsDataVector<MODEL, DATA>;
  template <typename DATA> using ObsDataPtr_ = boost::shared_ptr<ObsDataVector<MODEL, DATA> >;

 public:
// -----------------------------------------------------------------------------

  LETKF() {
    instantiateObsErrorFactory<MODEL>();
    instantiateObsFilterFactory<MODEL>();
  }

// -----------------------------------------------------------------------------

  virtual ~LETKF() {}

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const {
    // Setup observation window
    const eckit::LocalConfiguration windowConfig(fullConfig, "Assimilation Window");
    const util::Duration winlen(windowConfig.getString("Length"));
    const util::DateTime winbgn(windowConfig.getString("Begin"));
    const util::DateTime winend(winbgn + winlen);
    const util::DateTime winhalf = winbgn + winlen/2;
    Log::debug() << "Observation window is: " << windowConfig << std::endl;

    // Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "Geometry");
    const Geometry_ resol(resolConfig);

    // Setup model
    const eckit::LocalConfiguration modelConfig(fullConfig, "Model");
    const Model_ model(resol, modelConfig);

    // Setup observations
    const eckit::LocalConfiguration obsConfig(fullConfig, "Observations");
    Log::debug() << "Observation configuration is: " << obsConfig << std::endl;
    ObsSpaces_ obsdb(obsConfig, winbgn, winend);
    ObsAuxCtrls_ ybias(obsConfig);
    Observations_ yobs(obsdb, "ObsValue");

    // Setup initial obs error / qc
    // TODO(Travis) setup QC flags as well
    std::vector<ObsDataPtr_<float> > obserr;
    for (unsigned jj=0; jj < obsdb.size(); ++jj) {
      ObsDataPtr_<float> tmperr(new ObsData_<float>(
        obsdb[jj], obsdb[jj].obsvariables(), "ObsError"));
      obserr.push_back(tmperr);
      Log::info() << "initialize obs error: " << *tmperr << std::endl;
      obserr[jj]->save("EffectiveError");
    }
    ObsErrors_ rmat(obsConfig, obsdb);

    // Get initial state configurations
    const eckit::LocalConfiguration initialConfig(fullConfig, "Initial Condition");
    std::vector<eckit::LocalConfiguration> memberConfig;
    initialConfig.get("state", memberConfig);
    Log::debug() << "LETKF using " << memberConfig.size() << " members." << std::endl;

    // Loop over all ensemble members
    // TODO(Travis) use a StateEnsemble class
    ObsEnsemble_ obsens(obsdb, memberConfig.size());
    std::vector< std::shared_ptr<State_> > ens_xx;
    for (unsigned jj = 0; jj < memberConfig.size(); ++jj) {
      // Setup initial state
      // TODO(Travis) change the way input file name is specified, make
      //  more similar to how the output ens config is done
      Log::info() << std::endl << "Initial configuration for member " << jj+1 << " is: "
                  << memberConfig[jj] << std::endl;
      std::shared_ptr<State_> xx( new State_(resol, model.variables(), memberConfig[jj]));
      ModelAux_ moderr(resol, memberConfig[jj]);
      Log::test() << "Initial state for member " << jj+1 << ":" << *xx << std::endl;
      ens_xx.push_back(xx);

      // setup postprocessor: observers
      // TODO(Travis) obs filters not currently being used here, need to be added
      PostProcessor<State_> post;
      boost::shared_ptr<Observers<MODEL, State_> >
      pobs(new Observers<MODEL, State_> (obsConfig, obsdb, ybias));
      post.enrollProcessor(pobs);

      // compute H(x)
      // TODO(Travis) this is the 3D LETKF case, need to add model forecast for 4D LETKF
      post.initialize(*xx, winhalf, winlen);
      post.process(*xx);
      post.finalize(*xx);

      // save H(x)
      std::unique_ptr<Observations_> yeqv(pobs->release());
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << *yeqv << std::endl;
      yeqv->save("hofx0_"+std::to_string(jj+1));
      obsens[jj] = *yeqv;
    }

    // calculate background mean
    // TODO(Travis) make changes to IncrementEnsemble class so we can use that instead
    Accumulator<MODEL, State_, State_> bkg_mean(*ens_xx[0]);
    for (unsigned jj = 0; jj < ens_xx.size(); ++jj) {
      const double rr = 1.0/ens_xx.size();
      bkg_mean.accumul(rr, *ens_xx[jj]);
    }
    Log::test() << "Background mean :" << bkg_mean << std::endl;

    // calculate background ensemble perturbations
    // TODO(Travis) use the IncrementEnsemble class
    std::vector<std::shared_ptr<Increment_> > bkg_pert;
    for (unsigned jj = 0; jj < ens_xx.size(); ++jj) {
      std::shared_ptr<Increment_> dx (new Increment_(resol, model.variables(),
                                                     bkg_mean.validTime()));
      dx->diff(*ens_xx[jj], bkg_mean);
      bkg_pert.push_back(dx);
    }

    // TODO(Travis) optionally save the background mean / standard deviation

    // obs filters
    // TODO(Travis)

    // calculate H(x) ensemble mean
    // TODO(Travis) this should be pulled out and added as a method to ObsEnsemble
    Observations_ yb_mean(obsens[0]);
    for (unsigned jj = 1; jj < obsens.size(); ++jj) {
      Departures_ d(obsens[jj] - yb_mean);
      d *= 1.0/(jj+1.0);
      yb_mean += d;
    }
    Log::test() << "H(x) ensemble background mean: " << std::endl << yb_mean << std::endl;

    // calculate H(x) ensemble perturbations
    // TODO(Travis) pull this out into a DepartureEnsemble class
    std::vector< std::unique_ptr<Departures_> > ens_Yb;
    for (unsigned jj = 0; jj < obsens.size(); ++jj) {
      ens_Yb.push_back(std::unique_ptr<Departures_> (new Departures_(obsens[jj] - yb_mean)));
    }

    // calculate obs departures
    Departures_ ombg(yobs - yb_mean);
    ombg.save("ombg");
    Log::test() << "background y - H(x): " << std::endl << ombg << std::endl;

    // run the LETKF solver
    // TODO(Travis) put this inside a geometry iterator loop and run
    //  once for each gridpoint, with obs localization
    Log::info() << "Beginning core LETKF solver..." << std::endl;
    const eckit::LocalConfiguration letkfConfig(fullConfig, "letkf");
    const Eigen::MatrixXd trans = calcTrans(letkfConfig, ombg, ens_Yb, rmat);
    Log::info() << "LETKF solver completed." << std::endl;

    // use the transform matrix to calculate the analysis perturbations
    std::vector<std::shared_ptr<Increment_> > ana_pert;
    for (unsigned jj=0; jj < obsens.size(); ++jj) {
      std::shared_ptr<Increment_> dx(new Increment_(*bkg_pert[0]));
      (*dx) *= trans(0, jj);
      for (unsigned ii=1; ii < obsens.size(); ++ii) {
        Increment_ dx2 = *bkg_pert[ii];
        dx2 *= trans(ii, jj);
        (*dx) += dx2;
      }
      ana_pert.push_back(dx);
    }

    // calculate final analysis states
    for (unsigned jj = 0; jj < ens_xx.size(); ++jj) {
      *ens_xx[jj] = bkg_mean;
      *ens_xx[jj] += (*ana_pert[jj]);
    }

    // TODO(Travis) optionally save analysis mean / standard deviation

    // save the analysis ensemble
    for (unsigned jj=0; jj < obsens.size(); ++jj) {
      eckit::LocalConfiguration outConfig(fullConfig, "output");
      outConfig.set("member", jj+1.0);
      ens_xx[jj]->write(outConfig);
    }

    // calculate oman
    // TODO(Travis) redundant code with the previous H(x) loop... simplify this
    for (unsigned jj = 0; jj < memberConfig.size(); ++jj) {
      // setup postprocessor: observers
      PostProcessor<State_> post;
      boost::shared_ptr<Observers<MODEL, State_> >
      pobs(new Observers<MODEL, State_> (obsConfig, obsdb, ybias));
      post.enrollProcessor(pobs);

      // compute H(x)
      post.initialize(*ens_xx[jj], winhalf, winlen);
      post.process(*ens_xx[jj]);
      post.finalize(*ens_xx[jj]);

      // save H(x)
      std::unique_ptr<Observations_> yeqv(pobs->release());
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << *yeqv << std::endl;
      yeqv->save("hofx1_"+std::to_string(jj+1));
      obsens[jj] = *yeqv;
    }

    // calcualate H(x) ensemble analysis mean
    // TODO(Travis) this should be pulled out and added as a method to ObsEnsemble
    Observations_ ya_mean(obsens[0]);
    for (unsigned jj = 1; jj < obsens.size(); ++jj) {
      Departures_ d(obsens[jj] - ya_mean);
      d *= 1.0/(jj+1.0);
      ya_mean += d;
    }
    Log::test() << "H(x) ensemble analysis mean: " << std::endl << ya_mean << std::endl;

    // calculate analysis obs departures
    Departures_ oman(yobs - ya_mean);
    oman.save("oman");
    Log::test() << "analysis y - H(x): " << std::endl << oman << std::endl;

    // display overall background/analysis RMS stats
    double rms_bkg = sqrt(ombg.dot_product_with(ombg) / ombg.size());
    double rms_ana = sqrt(oman.dot_product_with(oman) / oman.size());
    Log::test() << "ombg RMS: " << rms_bkg << std::endl
                << "oman RMS: " << rms_ana << std::endl;

    return 0;
  }

// -----------------------------------------------------------------------------

 private:
  std::string appname() const {
    return "oops::LETKF<" + MODEL::name() + ">";
  }

// -----------------------------------------------------------------------------

/// Calculate the transform matrix to go from background to analysis perturbations
/*!
 * @param  conf  The LETKF portion of the yaml configuration file
 * @param  dy    Observation departures from background
 * @param  Yb    Background ensemble perturbations in observation space
 * @param  R     Observation error
 *
 * The following steps are performed to calculate the transform matix, \f$ T \f$
 *
 * - \f$ work = Y_b^T R^{-1} Y_b + (k-1)/\rho \f$
 * - eigenvector decomposition of \f$ work \f$. (used in the next two steps)
 * - \f$ \tilde{P_a}= [ work]^{-1}\f$
 * - \f$ W_a = [ (k-1) \tilde{P_a} ]^{1/2}\f$
 * - \f$ \bar{w_a} = \tilde{P_a} Y_b^T R^{-1} \delta y \f$
 * - \f$ \bar{w_a} \f$ is added to each column of \f$ W_a \f$ to form \f$ T \f$
 */
  static const Eigen::MatrixXd calcTrans(
                                  const eckit::Configuration &conf,
                                  const Departures_ & dy,
                                  const std::vector< std::unique_ptr<Departures_> > & Yb,
                                  const ObsErrors<MODEL> & R) {
    unsigned int nbv = Yb.size();  // number of ensemble members
    double infl = 1.0;    // TODO(Travis): read multiplicative inflation from config

    Eigen::MatrixXd work(nbv, nbv);
    Eigen::MatrixXd trans;

    // fill in the work matrix (note that since the matrix is symmetric,
    // only lower triangular half is filled)
    // work = Y^T R^-1 Y + (nbv-1)/infl I
    for (unsigned jj=0; jj < nbv; ++jj) {
      Departures_ Cj(*Yb[jj]);
      R.inverseMultiply(Cj);
      for (unsigned ii=jj; ii < nbv; ++ii) {
        work(ii, jj) = Cj.dot_product_with(*Yb[ii]);
      }
      work(jj, jj) += (nbv-1.0) / infl;
    }

    // eigenvalues and eigenvectors of the above matrix
    Eigen::VectorXd eival(nbv);
    Eigen::MatrixXd eivec(nbv, nbv);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(work);
    eival = es.eigenvalues().real();
    eivec = es.eigenvectors().real();

    // Pa   = [ Yb^T R^-1 Yb + (nbv-1)/infl I ] ^-1
    work = eivec * eival.cwiseInverse().asDiagonal() * eivec.transpose();

    // Wa = sqrt[ (nbv-1) Pa ]
    trans = eivec
      * ((nbv-1) * eival.array().inverse()).sqrt().matrix().asDiagonal()
      * eivec.transpose();

    // wa = Pa Yb^T R^-1 dy
    Eigen::VectorXd wa(nbv);
    Departures_ Rinvdy(dy);
    R.inverseMultiply(Rinvdy);
    for (unsigned jj=0; jj < nbv; ++jj) {
      wa(jj) = Yb[jj]->dot_product_with(Rinvdy);
    }
    wa = work * wa;

    // add wa to each column of Wa to get T
    trans.colwise() += wa;

    return trans;
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_LETKF_H_