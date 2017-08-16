/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_DIRAC_H_
#define OOPS_RUNS_DIRAC_H_

#include <sstream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "util/Logger.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/StateWriter.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/runs/Application.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

template <typename MODEL> class Dirac : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>      ModelAux_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;
  typedef Variables<MODEL>           Variables_;
  typedef Localization<MODEL>        Localization_;

 public:
// -----------------------------------------------------------------------------
  Dirac() {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~Dirac() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {

//  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);
    Log::info() << "Setup resolution OK" << std::endl;

//  Setup variables
    const eckit::LocalConfiguration varConfig(fullConfig, "variables");
    const Variables_ vars(varConfig);
    Log::info() << "Setup variables OK" << std::endl;

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "initial");
    const State_ xx(resol, initialConfig);
    Log::info() << "Setup initial state OK" << std::endl;

//  Setup times
    const util::DateTime bgndate(xx.validTime());
    Log::info() << "Setup times OK" << std::endl;

//  Setup localization
    const eckit::LocalConfiguration covarConfig(fullConfig, "Covariance");
    const eckit::LocalConfiguration locConfig(covarConfig, "localization");
    boost::scoped_ptr<Localization_> loc_;
    loc_.reset(new Localization_(xx, locConfig));
    Log::info() << "Setup localization OK" << std::endl;

//  Setup increment
    Increment_ dxinit(resol, vars, bgndate);
    const eckit::LocalConfiguration diracConfig(fullConfig, "dirac");
    dxinit.dirac(diracConfig);
    Log::info() << "Setup increment OK" << std::endl;

//  Apply NICAS
    Increment_ dx(dxinit);
    loc_->multiply(dx);
    Log::info() << "Apply NICAS OK" << std::endl;
    Log::test() << "Increment norm: " << dx.norm() << std::endl;

//  Write increment
    const eckit::LocalConfiguration output_nicas(fullConfig, "output_nicas");
    dx.write(output_nicas);
    Log::info() << "Write increment OK" << std::endl;

//  Setup full ensemble B matrix
    boost::scoped_ptr< ModelSpaceCovarianceBase<MODEL> >
      Bens(CovarianceFactory<MODEL>::create(covarConfig, resol, vars, xx));
    Bens->linearize(xx, resol);
    Log::info() << "Setup full ensemble B matrix OK" << std::endl;

//  Apply full ensemble B matrix
    Increment_ dxin(dxinit);
    Increment_ dxout(resol, vars, bgndate);
    Bens->multiply(dxin,dxout);
    Log::info() << "Apply full ensemble B matrix OK" << std::endl;
    Log::test() << "Increment norm: " << dxout.norm() << std::endl;

//  Write increment
    const eckit::LocalConfiguration output_Bens(fullConfig, "output_Bens");
    dxout.write(output_Bens);
    Log::info() << "Write increment OK" << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::Dirac<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_DIRAC_H_
