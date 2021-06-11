/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/L95Traits.h"
#include "oops/runs/Run.h"
#include "test/interface/ObsLocalization.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::ObsLocalization<lorenz95::L95Traits, lorenz95::L95ObsTraits> tests;
  return run.execute(tests);
}

