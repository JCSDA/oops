/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_EXTERNALDFI_H_
#define OOPS_RUNS_EXTERNALDFI_H_

#include <string>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/base/StateWriter.h"
#include "oops/base/WeightedMean.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Model.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/abor1_cpp.h"

namespace oops {

template <typename MODEL> class ExternalDFI : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef Model<MODEL>               Model_;
  typedef ModelAuxControl<MODEL>     ModelAux_;
  typedef State<MODEL>               State_;
 public:
// -----------------------------------------------------------------------------
  ExternalDFI() {}
// -----------------------------------------------------------------------------
  virtual ~ExternalDFI() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig);

//  Setup Model
    const eckit::LocalConfiguration modelConfig(fullConfig, "model");
    const Model_ model(resol, modelConfig);

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "initial");
    State_ xx(resol, initialConfig);
    Log::test() << "Initial state: " << xx.norm() << std::endl;

//  Setup augmented state
    const ModelAux_ moderr(resol, initialConfig);

//  Setup times
    const util::Duration fclength(fullConfig.getString("forecast_length"));
    const util::DateTime bgndate(xx.validTime());
    const util::DateTime enddate(bgndate + fclength);
    Log::info() << "Running forecast from " << bgndate << " to " << enddate << std::endl;

//  Setup post-processing
    PostProcessor<State_> post;

    eckit::LocalConfiguration prtConf;
    fullConfig.get("prints", prtConf);
    post.enrollProcessor(new StateInfo<State_>("fc", prtConf));

//  Setup DFI
    PostProcessor<State_> pp(post);

    const eckit::LocalConfiguration dfiConf(fullConfig, "dfi");
    const util::Duration dfispan(dfiConf.getString("filter_span"));
    const util::DateTime dfitime(bgndate+dfispan/2);
    boost::shared_ptr< WeightedMean<MODEL, State_> >
      pdfi(new WeightedMean<MODEL, State_>(dfitime, dfispan, resol, dfiConf));
    pp.enrollProcessor(pdfi);

//  Run DFI forecast
    model.forecast(xx, moderr, dfispan, pp);

//  Retrieve initialized state
    boost::scoped_ptr<State_> xdfi(pdfi->releaseMean());
    Log::test() << "Filtered state: " << xdfi->norm() << std::endl;

//  Setup forecast outputs
    const eckit::LocalConfiguration outConfig(fullConfig, "output");
    post.enrollProcessor(new StateWriter<State_>(bgndate, outConfig));

//  Run forecast from initialized state
    const util::Duration fclen = fclength - dfispan/2;
    if (fclength < util::Duration(0)) ABORT("DFI: fliter span longer than forecast");
    State_ zz(*xdfi);
    model.forecast(zz, moderr, fclen, post);
    Log::test() << "Final state: " << zz.norm() << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::ExternalDFI<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_EXTERNALDFI_H_
