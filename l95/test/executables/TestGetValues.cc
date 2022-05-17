/*
 * (C) Copyright 2020-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/L95Traits.h"
#include "oops/runs/Run.h"
#include "test/interface/GetValues.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::GetValues<lorenz95::L95Traits, lorenz95::L95ObsTraits> tests;
  return run.execute(tests);
}
