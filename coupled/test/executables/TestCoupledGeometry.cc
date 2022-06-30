/*
 * (C) Copyright 2021- UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "lorenz95/L95Traits.h"
#include "model/QgTraits.h"
#include "oops/coupled/TraitCoupled.h"
#include "oops/runs/Run.h"
#include "test/interface/Geometry.h"

int main(int argc, char ** argv) {
  oops::Run run(argc, argv);
  test::Geometry<oops::TraitCoupled<qg::QgTraits, lorenz95::L95Traits>> tests;
  return run.execute(tests);
}
