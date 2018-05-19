/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_BASE_INSTANTIATEFILTERFACTORY_H_
#define OOPS_BASE_INSTANTIATEFILTERFACTORY_H_

#include "oops/base/FilterBase.h"
#include "oops/base/GeoVaLsWriter.h"

namespace oops {

template <typename MODEL> void instantiateFilterFactory() {
  static FilterMaker<MODEL, GeoVaLsWriter<MODEL> >   makerGVWriter_("GOMsaver");
}

}  // namespace oops

#endif  // OOPS_BASE_INSTANTIATEFILTERFACTORY_H_
