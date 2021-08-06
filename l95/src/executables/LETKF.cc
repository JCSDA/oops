/*
 * (C) Copyright 2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/instantiateL95ChangeVarFactory.h"
#include "lorenz95/L95Traits.h"
#include "oops/runs/LocalEnsembleDA.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  lorenz95::instantiateL95ChangeVarFactory();
  oops::LocalEnsembleDA<lorenz95::L95Traits, lorenz95::L95ObsTraits> letkf;
  return run.execute(letkf);
}
