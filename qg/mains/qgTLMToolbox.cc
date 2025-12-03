/*
 * (C) Copyright 2023-2025 UCAR.
 * (C) Crown copyright 2023-2025 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/QgTraits.h"

#include "oops/runs/Run.h"
#include "oops/runs/TLMToolbox.h"

int main(int argc, char ** argv) {
  oops::Run run(argc, argv);
  oops::TLMToolbox<qg::QgTraits> tlmToolbox;
  return run.execute(tlmToolbox);
}
