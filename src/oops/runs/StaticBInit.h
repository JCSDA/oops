/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_STATICBINIT_H_
#define OOPS_RUNS_STATICBINIT_H_

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class StaticBInit : public Application {
  typedef ModelSpaceCovarianceBase<MODEL>  Covariance_;
  typedef Geometry<MODEL>                  Geometry_;
  typedef Increment<MODEL>                 Increment_;
  typedef State<MODEL>                     State_;

 public:
  // -----------------------------------------------------------------------------
  explicit StaticBInit(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }
  // -----------------------------------------------------------------------------
  virtual ~StaticBInit() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    //  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "geometry");
    const Geometry_ resol(resolConfig, this->getComm());

    //  Setup variables
    const Variables vars(fullConfig, "analysis variables");

    //  Setup background state
    const eckit::LocalConfiguration bkgconf(fullConfig, "background");
    State_ xx(resol, bkgconf);

    //  Initialize static B matrix
    const eckit::LocalConfiguration covarconf(fullConfig, "background error");
    std::unique_ptr< Covariance_ >
     Bmat(CovarianceFactory<MODEL>::create(covarconf, resol, vars, xx, xx));

    //  Randomize B matrix
    Increment_ dx(resol, vars, xx.validTime());
    Bmat->randomize(dx);
    Log::test() << dx << std::endl;

    return 0;
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::StaticBInit<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_STATICBINIT_H_
