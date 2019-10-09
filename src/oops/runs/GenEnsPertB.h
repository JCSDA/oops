/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_GENENSPERTB_H_
#define OOPS_RUNS_GENENSPERTB_H_

#include <memory>
#include <sstream>
#include <string>


#include "eckit/config/Configuration.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateWriter.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class GenEnsPertB : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>      ModelAux_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
// -----------------------------------------------------------------------------
  explicit GenEnsPertB(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~GenEnsPertB() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig, this->getComm());

//  Setup variables
    const Variables vars(fullConfig);

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "initial");
    const State_ xx(resol, model.variables(), initialConfig);
    Log::test() << "Initial state: " << xx << std::endl;

//  Setup augmented state
    const ModelAux_ moderr(resol, initialConfig);

//  Setup times
    const util::Duration fclength(fullConfig.getString("forecast_length"));
    const util::DateTime bgndate(xx.validTime());
    const util::DateTime enddate(bgndate + fclength);
    Log::info() << "Running forecast from " << bgndate << " to " << enddate << std::endl;

//  Setup B matrix
    const eckit::LocalConfiguration covar(fullConfig, "Covariance");
    std::unique_ptr< ModelSpaceCovarianceBase<MODEL> >
      Bmat(CovarianceFactory<MODEL>::create(covar, resol, vars, xx, xx));

//  Generate perturbed states
    Increment_ dx(resol, vars, bgndate);
    const int members = fullConfig.getInt("members");
    for (int jm = 0; jm < members; ++jm) {
//    Generate pertubation
      Bmat->randomize(dx);

//    Add mean state
      State_ xp(xx);
      xp += dx;

//    Setup forecast outputs
      PostProcessor<State_> post;

      eckit::LocalConfiguration outConfig(fullConfig, "output");
      outConfig.set("member", jm+1);

      post.enrollProcessor(new StateWriter<State_>(bgndate, outConfig));

//    Run forecast
      model.forecast(xp, moderr, fclength, post);
      Log::test() << "Member " << jm << " final state: " << xp << std::endl;
    }

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::GenEnsPertB<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_GENENSPERTB_H_
