/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_INSTANTIATECOSTFACTORY_H_
#define OOPS_ASSIMILATION_INSTANTIATECOSTFACTORY_H_

#include "oops/assimilation/CostFct3DVar.h"
#include "oops/assimilation/CostFct4DEnsVar.h"
#include "oops/assimilation/CostFct4DVar.h"
#include "oops/assimilation/CostFctFGAT.h"
#include "oops/assimilation/CostFctWeak.h"
#include "oops/assimilation/CostFunction.h"

namespace oops {

template <typename MODEL, typename OBS> void instantiateCostFactory() {
  static CostMaker<MODEL, OBS, CostFct3DVar<MODEL, OBS> >    maker3DVar_("3D-Var");
  static CostMaker<MODEL, OBS, CostFct4DVar<MODEL, OBS> >    maker4DVar_("4D-Var");
  static CostMaker<MODEL, OBS, CostFct4DEnsVar<MODEL, OBS> > maker4DEns_("4D-Ens-Var");
  static CostMaker<MODEL, OBS, CostFctFGAT<MODEL, OBS> >     maker3FGAT_("3D-FGAT");
  static CostMaker<MODEL, OBS, CostFctWeak<MODEL, OBS> >     maker4Weak_("4D-Weak");
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_INSTANTIATECOSTFACTORY_H_
