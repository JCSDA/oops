/*
 * (C) Copyright 2025- UCAR
 * (C) Copyright 2025- NOAA/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/QgObsTraits.h"
#include "model/QgTraits.h"

#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/coupled/BlockDiagonalCovarianceCoupled.h"
#include "oops/coupled/GetValuesCoupled.h"
#include "oops/coupled/TraitCoupled.h"

#include "oops/runs/Run.h"
#include "oops/runs/Variational.h"

#include "./QgTraits2.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::instantiateCovarFactory<qg::QgTraits>();
  oops::instantiateCovarFactory<qg::QgTraits2>();
  oops::instantiateCovarFactory<oops::TraitCoupled<qg::QgTraits, qg::QgTraits2> >();
  static oops::CovarMaker<oops::TraitCoupled<qg::QgTraits, qg::QgTraits2>,
      oops::BlockDiagonalCovarianceCoupled<qg::QgTraits, qg::QgTraits2> >
        makerCoupled_("Coupled Block Diagonal");
  oops::Variational<oops::TraitCoupled<qg::QgTraits, qg::QgTraits2>, qg::QgObsTraits> var;
  return run.execute(var);
}
