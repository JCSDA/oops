/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/QgTraits.h"
#include "oops/runs/Run.h"
#include "test/base/ModelIncrement.h"

int main(const int argc, const char ** argv) {
  oops::Run run(argc, argv);
  test::ModelIncrement<qg::QgTraits> tests;
  run.execute(tests);
  return 0;
};

