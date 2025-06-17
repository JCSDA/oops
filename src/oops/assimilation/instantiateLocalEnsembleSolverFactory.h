/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_INSTANTIATELOCALENSEMBLESOLVERFACTORY_H_
#define OOPS_ASSIMILATION_INSTANTIATELOCALENSEMBLESOLVERFACTORY_H_

#include "oops/assimilation/GETKFSolver.h"
#include "oops/assimilation/GETKFSolverPert.h"
#include "oops/assimilation/LETKFSolver.h"
#include "oops/assimilation/LETKFSolverPert.h"
#include "oops/assimilation/LocalEnsembleSolver.h"

namespace oops {

template <typename MODEL, typename OBS> void instantiateLocalEnsembleSolverFactory() {
  static LocalEnsembleSolverMaker<MODEL, OBS, DeterministicLETKF<MODEL, OBS> >
    makerLETKF_("Deterministic LETKF");
  static LocalEnsembleSolverMaker<MODEL, OBS, StochasticLETKF<MODEL, OBS> >
    makerLETKFPert_("Stochastic LETKF");
  static LocalEnsembleSolverMaker<MODEL, OBS, DeterministicGETKF<MODEL, OBS> >
    makerGETKF_("Deterministic GETKF");
  static LocalEnsembleSolverMaker<MODEL, OBS, StochasticGETKF<MODEL, OBS> >
    makerGETKFPert_("Stochastic GETKF");
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_INSTANTIATELOCALENSEMBLESOLVERFACTORY_H_
