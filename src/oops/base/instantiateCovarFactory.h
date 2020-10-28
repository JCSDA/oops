/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_INSTANTIATECOVARFACTORY_H_
#define OOPS_BASE_INSTANTIATECOVARFACTORY_H_

#include "oops/base/EnsembleCovariance.h"
#include "oops/base/HybridCovariance.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/interface/ErrorCovariance.h"

namespace oops {

template <typename MODEL> void instantiateCovarFactory() {
  static CovarMaker<MODEL, EnsembleCovariance<MODEL> >  makerEnsemble_("ensemble");
  static CovarMaker<MODEL, HybridCovariance<MODEL> >    makerHybrid_("hybrid");
  static CovarMaker<MODEL, ErrorCovariance<MODEL> >     makerModel_(MODEL::nameCovar());
}

}  // namespace oops

#endif  // OOPS_BASE_INSTANTIATECOVARFACTORY_H_
