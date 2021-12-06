/*
 * (C) Copyright 2021- UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "lorenz95/L95Traits.h"
#include "model/QgTraits.h"

#include "oops/coupled/instantiateCoupledFactory.h"
#include "oops/coupled/TraitCoupled.h"
#include "oops/generic/instantiateModelFactory.h"
#include "oops/runs/Run.h"

#include "test/interface/Model.h"

int main(int argc, char ** argv) {
  oops::Run run(argc, argv);
  oops::instantiateModelFactory<qg::QgTraits>();
  oops::instantiateModelFactory<lorenz95::L95Traits>();
  oops::instantiateCoupledFactory<qg::QgTraits, lorenz95::L95Traits>();
  test::Model<oops::TraitCoupled<qg::QgTraits, lorenz95::L95Traits>> tests;
  return run.execute(tests);
}
