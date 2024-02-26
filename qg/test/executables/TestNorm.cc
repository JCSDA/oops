/*
 * (C) Crown Copyright 2023, the Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/instantiateQgNormFactory.h"
#include "model/QgTraits.h"
#include "oops/runs/Run.h"
#include "test/base/NormBase.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  qg::instantiateQgNormFactory();
  test::TestNorm<qg::QgTraits> tests;
  return run.execute(tests);
}
