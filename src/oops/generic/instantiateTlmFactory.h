/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_GENERIC_INSTANTIATETLMFACTORY_H_
#define OOPS_GENERIC_INSTANTIATETLMFACTORY_H_

#include "oops/generic/LinearModelId.h"
#include "oops/interface/LinearModelBase.h"

namespace oops {

template <typename MODEL> void instantiateTlmFactory() {
  static LinearModelMaker<MODEL, LinearModelId<MODEL> >  makerTLMID_("identity");
}

}  // namespace oops

#endif  // OOPS_GENERIC_INSTANTIATETLMFACTORY_H_
