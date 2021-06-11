/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/instantiateL95ChangeVarFactory.h"
#include "lorenz95/L95Traits.h"
#include "oops/runs/HofX3D.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  lorenz95::instantiateL95ChangeVarFactory();
  oops::HofX3D<lorenz95::L95Traits, lorenz95::L95ObsTraits> hofx;
  return run.execute(hofx);
}
