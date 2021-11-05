/*
 * (C) Copyright 2009-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include "lorenz95/L95Traits.h"
#include "oops/runs/Run.h"
#include "oops/runs/SqrtOfVertLoc.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::SqrtOfVertLoc<lorenz95::L95Traits> sqrtgen;
  return run.execute(sqrtgen);
}
