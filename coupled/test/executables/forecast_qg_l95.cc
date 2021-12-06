/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/L95Traits.h"
#include "oops/qg/QgTraits.h"

#include "oops/generic/instantiateModelFactory.h"
#include "oops/runs/Forecast.h"
#include "oops/runs/Run.h"

#include "oops/coupled/instantiateCoupledFactory.h"
#include "oops/coupled/TraitCoupled.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);

  oops::instantiateModelFactory<qg::QgTraits>();
  oops::instantiateModelFactory<lorenz95::L95Traits>();
  oops::instantiateCoupledFactory<qg::QgTraits, lorenz95::L95Traits>();

  oops::Forecast<oops::TraitCoupled<qg::QgTraits, lorenz95::L95Traits>> fc;

  return run.execute(fc);
}

