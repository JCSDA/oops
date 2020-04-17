/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/instantiateQgChangeVarFactory.h"
#include "model/instantiateQgLocalizationFactory.h"
#include "model/QgTraits.h"
#include "oops/runs/EDA.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  qg::instantiateQgChangeVarFactory();
  qg::instantiateQgLocalizationFactory();
  oops::EDA<qg::QgTraits> eda;
  return run.execute(eda);
}
