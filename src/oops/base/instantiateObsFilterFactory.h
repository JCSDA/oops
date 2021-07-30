/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_INSTANTIATEOBSFILTERFACTORY_H_
#define OOPS_BASE_INSTANTIATEOBSFILTERFACTORY_H_

#include "oops/generic/GeoVaLsWriter.h"
#include "oops/generic/ObsFilterBase.h"

namespace oops {

template <typename OBS> void instantiateObsFilterFactory() {
  static FilterMaker<OBS, GeoVaLsWriter<OBS>> makerGVWriter_("GOMsaver");
}

}  // namespace oops

#endif  // OOPS_BASE_INSTANTIATEOBSFILTERFACTORY_H_
