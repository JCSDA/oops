/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_HOFX3D_H_
#define OOPS_RUNS_HOFX3D_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObsFilters.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/StateInfo.h"
#include "oops/generic/instantiateObsErrorFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class HofX3D : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef ObsAuxControls<MODEL>      ObsAuxCtrls_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsErrors<MODEL>           ObsErrors_;
  typedef ObsFilters<MODEL>          ObsFilters_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
  typedef State<MODEL>               State_;
  typedef boost::shared_ptr<ObsFilters_> PtrFilters_;
  template <typename DATA> using ObsData_ = ObsDataVector<MODEL, DATA>;
  template <typename DATA> using ObsDataPtr_ = boost::shared_ptr<ObsData_<DATA> >;

 public:
// -----------------------------------------------------------------------------
  HofX3D() {
    instantiateObsFilterFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofX3D() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const eckit::LocalConfiguration windowConf(fullConfig, "Assimilation Window");
    const util::Duration winlen(windowConf.getString("Length"));
    const util::DateTime winbgn(windowConf.getString("Begin"));
    const util::DateTime winend(winbgn + winlen);
    const util::DateTime winhalf = winbgn + winlen/2;
    Log::info() << "Observation window is:" << windowConf << std::endl;

//  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "Geometry");
    const Geometry_ resol(resolConfig);

//  Setup initial state
    const eckit::LocalConfiguration initialConfig(fullConfig, "Initial Condition");
    Log::info() << "Initial configuration is:" << initialConfig << std::endl;
    oops::Variables vars(initialConfig);
    State_ xx(resol, vars, initialConfig);
    Log::test() << "Initial state: " << xx << std::endl;

//  Setup forecast outputs
    PostProcessor<State_> post;

    eckit::LocalConfiguration prtConf;
    fullConfig.get("Prints", prtConf);
    post.enrollProcessor(new StateInfo<State_>("fc", prtConf));

//  Setup observations
    const eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    Log::info() << "Observations configuration is:" << obsconf << std::endl;
    ObsSpaces_ obspace(obsconf, winbgn, winend);

//  Setup observations bias
    ObsAuxCtrls_ ybias(obsconf);

//  Setup QC filters
    std::vector<eckit::LocalConfiguration> typeconfs;
    obsconf.get("ObsTypes", typeconfs);
    std::vector<ObsDataPtr_<int> > qcflags_;
    std::vector<ObsDataPtr_<float> > obserr_;
    std::vector<PtrFilters_> filters;

    for (size_t jj = 0; jj < obspace.size(); ++jj) {
//    Allocate QC flags
      ObsDataPtr_<int> tmpqc(new ObsData_<int>(obspace[jj], obspace[jj].obsvariables()));
      qcflags_.push_back(tmpqc);

//    Allocate and read initial obs error
      ObsDataPtr_<float> tmperr(new ObsData_<float>(obspace[jj],
                                obspace[jj].obsvariables(), "ObsError"));
      obserr_.push_back(tmperr);

      PtrFilters_ tmp(new ObsFilters_(obspace[jj], typeconfs[jj], qcflags_[jj], obserr_[jj]));
      filters.push_back(tmp);
    }

//  Setup Observers
    boost::shared_ptr<Observers<MODEL, State_> >
      pobs(new Observers<MODEL, State_>(obsconf, obspace, ybias, filters));
    post.enrollProcessor(pobs);

//  Compute H(x)
    post.initialize(xx, winhalf, winlen);
    post.process(xx);
    post.finalize(xx);

//  Save QC flags
    for (size_t jj = 0; jj < obspace.size(); ++jj) {
      qcflags_[jj]->save("EffectiveQC");
      obserr_[jj]->save("EffectiveError");
    }

//  Save H(x)
    std::unique_ptr<Observations_> yobs(pobs->release());
    Log::test() << "H(x): " << *yobs << std::endl;
    yobs->save("hofx");

    post.clear();

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::HofX3D<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HOFX3D_H_
