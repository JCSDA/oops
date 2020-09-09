/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_RTPP_H_
#define OOPS_RUNS_RTPP_H_

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/StateEnsemble.h"
#include "oops/interface/Geometry.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"


namespace oops {

/// \brief Application for relaxation to prior perturbation (RTPP) inflation
template <typename MODEL> class RTPP : public Application {
  typedef Geometry<MODEL>                  Geometry_;
  typedef Increment4D<MODEL>               Increment4D_;
  typedef State4D<MODEL>                   State4D_;
  typedef StateEnsemble<MODEL>             StateEnsemble_;

 public:
// -----------------------------------------------------------------------------

  explicit RTPP(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}

// -----------------------------------------------------------------------------

  virtual ~RTPP() = default;

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const {
    // Setup geometry
    const eckit::LocalConfiguration geometryConfig(fullConfig, "geometry");
    const Geometry_ geometry(geometryConfig, this->getComm());

    // Get configurations
    const eckit::LocalConfiguration bgConfig(fullConfig, "background");
    const eckit::LocalConfiguration anConfig(fullConfig, "analysis");
    const float factor = fullConfig.getFloat("factor");

    // Read all ensemble members
    StateEnsemble_ bgens(geometry, bgConfig);
    StateEnsemble_ anens(geometry, anConfig);
    const size_t nens = bgens.size();
    ASSERT(nens == anens.size());

    Variables anvars = anens.variables();
    if (fullConfig.has("analysis variables"))
      anvars = Variables(fullConfig, "analysis variables");

    // calculate ensemble means
    State4D_ bg_mean = bgens.mean();
    State4D_ an_mean = anens.mean();

    Log::test() << "Background member 1:" << bgens[0] << std::endl;
    Log::test() << "Analysis member 1:" << anens[0] << std::endl;

    // update analysis with RTPP
    for (size_t jj = 0; jj < nens; ++jj) {
      // calculate RTPP perturbation
      Increment4D_ pertTot(geometry, anvars, anens[jj].validTimes());
      pertTot.zero();

      Increment4D_ pert(pertTot);

      pert.diff(bgens[jj], bg_mean);
      pertTot.axpy(factor, pert);

      pert.diff(anens[jj], an_mean);
      pertTot.axpy((1.0 - factor), pert);

      // an_mean contains copies of anens[0] for non-state variables
      // using zero+accumul instead of "=" ensures that only
      // state variables are modified in anens[jj]
      anens[jj].zero();
      anens[jj].accumul(1.0, an_mean);

      // add analysis variable perturbations
      anens[jj] += pertTot;
    }
    Log::test() << "Updated Analysis member 1:" << anens[0] << std::endl;

    // save the analysis mean
    an_mean = anens.mean();   // calculate analysis mean
    Log::test() << "Analysis mean:" << an_mean << std::endl;
    eckit::LocalConfiguration outConfig(fullConfig, "output");
    outConfig.set("member", 0);
    an_mean.write(outConfig);

    // save the analysis ensemble
    size_t mymember;
    for (size_t jj=0; jj < nens; ++jj) {
      mymember = jj+1;
      outConfig.set("member", mymember);
      anens[jj].write(outConfig);
    }

    return 0;
  }

// -----------------------------------------------------------------------------

 private:
  std::string appname() const {
    return "oops::RTPP<" + MODEL::name() + ">";
  }

// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_RTPP_H_
