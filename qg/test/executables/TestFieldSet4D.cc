/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/QgTraits.h"
#include "oops/runs/Run.h"
#include "test/interface/FieldSet4D.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::FieldSet4D<qg::QgTraits> tests;
  return run.execute(tests);
}
