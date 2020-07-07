/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_HOFXNOMODEL_H_
#define OOPS_RUNS_HOFXNOMODEL_H_

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CalcHofX.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/Geometry.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL, typename OBS> class HofXNoModel : public Application {
  typedef Geometry<MODEL>            Geometry_;
  typedef Observations<OBS>          Observations_;
  typedef ObsSpaces<OBS>             ObsSpaces_;
  typedef State4D<MODEL>             State4D_;

 public:
// -----------------------------------------------------------------------------
  explicit HofXNoModel(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateObsFilterFactory<OBS>();
  }
// -----------------------------------------------------------------------------
  virtual ~HofXNoModel() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const eckit::LocalConfiguration windowConf(fullConfig, "Assimilation Window");
    const util::Duration winlen(windowConf.getString("window_length"));
    const util::DateTime winbgn(windowConf.getString("window_begin"));
    const util::DateTime winend(winbgn + winlen);

//  Setup geometry
    const eckit::LocalConfiguration geometryConfig(fullConfig, "Geometry");
    const Geometry_ geometry(geometryConfig, this->getComm());

//  Setup states for H(x)
    const eckit::LocalConfiguration stateConfig(fullConfig, "Forecasts");
    Log::info() << "States configuration is:" << stateConfig << std::endl;
    oops::Variables vars(stateConfig);
    State4D_ xx(geometry, vars, stateConfig);
    Log::test() << "Initial state: " << xx[0] << std::endl;

//  Setup observations
    const eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    Log::info() << "Observations configuration is:" << obsconf << std::endl;
    ObsSpaces_ obspace(obsconf, this->getComm(), winbgn, winend);

//  Setup and run observer
    CalcHofX<MODEL, OBS> hofx(obspace, geometry, fullConfig);
    const Observations_ & yobs = hofx.compute(xx);
    for (size_t jj = 0; jj < obspace.size(); ++jj) {
      hofx.qcFlags(jj).save("EffectiveQC");
      hofx.obsErrors(jj).save("EffectiveError");
    }
    Log::test() << "Final state: " << xx[xx.size()-1] << std::endl;

//  Save H(x)
    Log::test() << "H(x): " << std::endl << yobs << "End H(x)" << std::endl;
    yobs.save("hofx");

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::HofXNoModel<" + MODEL::name() + ", " + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_HOFXNOMODEL_H_
