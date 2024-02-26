/*
 * (C) Crown Copyright 2023, the Met Office.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/qg/instantiateQgLocalizationFactory.h"
#include "oops/qg/QgTraits.h"
#include "oops/runs/ControlPert.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  qg::instantiateQgLocalizationFactory();
  oops::ControlPert<qg::QgTraits, qg::QgObsTraits> var;
  return run.execute(var);
}
