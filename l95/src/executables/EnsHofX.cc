/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/L95Traits.h"
#include "oops/runs/EnsembleApplication.h"
#include "oops/runs/HofX4D.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::EnsembleApplication<oops::HofX4D <lorenz95::L95Traits, lorenz95::L95ObsTraits> > enshofx;
  return run.execute(enshofx);
}
