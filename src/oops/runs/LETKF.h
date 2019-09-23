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
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/GridPoint.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
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
  typedef ModelAuxControl<MODEL>    ModelAux_;
  typedef Departures<MODEL>         Departures_;
  typedef DeparturesEnsemble<MODEL> DeparturesEnsemble_;
  typedef Geometry<MODEL>           Geometry_;
  typedef GeometryIterator<MODEL>   GeometryIterator_;
  typedef Increment<MODEL>          Increment_;
  typedef Model<MODEL>              Model_;
  typedef ObsAuxControls<MODEL>     ObsAuxCtrls_;
  typedef ObsEnsemble<MODEL>        ObsEnsemble_;
  typedef ObsErrors<MODEL>          ObsErrors_;
  typedef ObsSpaces<MODEL>          ObsSpaces_;
  typedef Observations<MODEL>       Observations_;
  typedef State<MODEL>              State_;
  typedef StateEnsemble<MODEL>      StateEnsemble_;
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

    // Loop over all ensemble members
    StateEnsemble_ ens_xx(resol, model.variables(), initialConfig);
    ObsEnsemble_ obsens(obsdb, ens_xx.size());
    for (unsigned jj = 0; jj < ens_xx.size(); ++jj) {
      // TODO(Travis) change the way input file name is specified, make
      //  more similar to how the output ens config is done
      Log::test() << "Initial state for member " << jj+1 << ":" << ens_xx[jj] << std::endl;

      // setup postprocessor: observers
      // TODO(Travis) obs filters not currently being used here, need to be added
      PostProcessor<State_> post;
      boost::shared_ptr<Observers<MODEL, State_> >
      pobs(new Observers<MODEL, State_> (obsConfig, obsdb, ybias));
      post.enrollProcessor(pobs);

      // compute H(x)
      // TODO(Travis) this is the 3D LETKF case, need to add model forecast for 4D LETKF
      post.initialize(ens_xx[jj], winhalf, winlen);
      post.process(ens_xx[jj]);
      post.finalize(ens_xx[jj]);

      // save H(x)
      std::unique_ptr<Observations_> yeqv(pobs->release());
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << *yeqv << std::endl;
      yeqv->save("hofx0_"+std::to_string(jj+1));
      obsens[jj] = *yeqv;
    }

    // calculate background mean
    State_ bkg_mean(ens_xx.mean());
    Log::test() << "Background mean :" << bkg_mean << std::endl;

    // calculate background ensemble perturbations
    // TODO(Travis) use the IncrementEnsemble class
    std::vector<std::shared_ptr<Increment_> > bkg_pert;
    for (unsigned jj = 0; jj < ens_xx.size(); ++jj) {
      std::shared_ptr<Increment_> dx (new Increment_(resol, model.variables(),
                                                     bkg_mean.validTime()));
      dx->diff(ens_xx[jj], bkg_mean);
      bkg_pert.push_back(dx);
    }

    // TODO(Travis) optionally save the background mean / standard deviation

    // obs filters
    // TODO(Travis)

    // calculate H(x) ensemble mean
    Observations_ yb_mean(obsens.mean());
    Log::test() << "H(x) ensemble background mean: " << std::endl << yb_mean << std::endl;

    // calculate H(x) ensemble perturbations
    DeparturesEnsemble_ ens_Yb(obsens, yb_mean);

    // calculate obs departures
    Departures_ ombg(yobs - yb_mean);
    ombg.save("ombg");
    Log::test() << "background y - H(x): " << std::endl << ombg << std::endl;

    // initialize empty analysis perturbations
    std::vector<std::shared_ptr<Increment_> > ana_pert;
    for (unsigned jj=0; jj < obsens.size(); ++jj) {
      std::shared_ptr<Increment_> dx(new Increment_(*bkg_pert[0], false));
      ana_pert.push_back(dx);
    }

    // run the LETKF solver at each gridpoint
    Log::info() << "Beginning core LETKF solver..." << std::endl;
    const eckit::LocalConfiguration letkfConfig(fullConfig, "letkf");
    for (GeometryIterator_ i = resol.begin(); i != resol.end(); ++i) {
      // TODO(Travis) Next step is to use localized observations
      const Eigen::MatrixXd trans = calcTrans(letkfConfig, ombg, ens_Yb, rmat);

      // use the transform matrix to calculate the analysis perturbations
      for (unsigned jj=0; jj < obsens.size(); ++jj) {
        GridPoint gp = bkg_pert[0]->getPoint(i);
        gp *= trans(0, jj);
        for (unsigned ii=1; ii < obsens.size(); ++ii) {
          GridPoint gp2 = bkg_pert[ii]->getPoint(i);
          gp2 *= trans(ii, jj);
          gp += gp2;
        }
        ana_pert[jj]->setPoint(gp, i);
      }
    }
    Log::info() << "LETKF solver completed." << std::endl;

    // calculate final analysis states
    for (unsigned jj = 0; jj < ens_xx.size(); ++jj) {
      ens_xx[jj] = bkg_mean;
      ens_xx[jj] += (*ana_pert[jj]);
    }

    // TODO(Travis) optionally save analysis mean / standard deviation

    // save the analysis ensemble
    for (unsigned jj=0; jj < obsens.size(); ++jj) {
      eckit::LocalConfiguration outConfig(fullConfig, "output");
      outConfig.set("member", jj+1.0);
      ens_xx[jj].write(outConfig);
    }

    // calculate oman
    // TODO(Travis) redundant code with the previous H(x) loop... simplify this
    for (unsigned jj = 0; jj < ens_xx.size(); ++jj) {
      // setup postprocessor: observers
      PostProcessor<State_> post;
      boost::shared_ptr<Observers<MODEL, State_> >
      pobs(new Observers<MODEL, State_> (obsConfig, obsdb, ybias));
      post.enrollProcessor(pobs);

      // compute H(x)
      post.initialize(ens_xx[jj], winhalf, winlen);
      post.process(ens_xx[jj]);
      post.finalize(ens_xx[jj]);

      // save H(x)
      std::unique_ptr<Observations_> yeqv(pobs->release());
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << *yeqv << std::endl;
      yeqv->save("hofx1_"+std::to_string(jj+1));
      obsens[jj] = *yeqv;
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
                                  const DeparturesEnsemble_ & Yb,
                                  const ObsErrors<MODEL> & R) {
    unsigned int nbv = Yb.size();  // number of ensemble members
    double infl = 1.0;    // TODO(Travis): read multiplicative inflation from config

    Eigen::MatrixXd work(nbv, nbv);
    Eigen::MatrixXd trans;

    // fill in the work matrix (note that since the matrix is symmetric,
    // only lower triangular half is filled)
    // work = Y^T R^-1 Y + (nbv-1)/infl I
    for (unsigned jj=0; jj < nbv; ++jj) {
      Departures_ Cj(Yb[jj]);
      R.inverseMultiply(Cj);
      for (unsigned ii=jj; ii < nbv; ++ii) {
        work(ii, jj) = Cj.dot_product_with(Yb[ii]);
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
      wa(jj) = Yb[jj].dot_product_with(Rinvdy);
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
