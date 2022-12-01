/*
 * (C) Copyright 2022 MetOffice.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/QgTraits.h"
#include "oops/runs/HTLMCoefficients.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::HTLMCoefficients<qg::QgTraits> coeffs;
  return run.execute(coeffs);
}

