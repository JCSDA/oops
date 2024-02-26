/*
 * (C) Copyright 2022-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/qg/QgTraits.h"

#include "oops/coupled/GetValuesCoupled.h"
#include "oops/coupled/TraitCoupled.h"

#include "oops/runs/HofX3D.h"
#include "oops/runs/Run.h"

#include "./QgTraits2.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::HofX3D<oops::TraitCoupled<qg::QgTraits, qg::QgTraits2>, qg::QgObsTraits> hofx;

  return run.execute(hofx);
}
