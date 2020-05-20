/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_INSTANTIATELETKFSOLVERFACTORY_H_
#define OOPS_ASSIMILATION_INSTANTIATELETKFSOLVERFACTORY_H_

#include "oops/assimilation/LETKFSolverBase.h"
#include "oops/assimilation/LETKFSolverOOPS.h"

namespace oops {

template <typename MODEL> void instantiateLETKFSolverFactory() {
//  static LETKFSolverMaker<MODEL, LETKFSolver<MODEL> > maker2_("LETKF");
  static LETKFSolverMaker<MODEL, LETKFSolverOOPS<MODEL> > maker1_("LETKF_OOPS");
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_INSTANTIATELETKFSOLVERFACTORY_H_
