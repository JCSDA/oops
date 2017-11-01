/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_INSTANTIATEQGOBSFACTORY_H_
#define QG_MODEL_INSTANTIATEQGOBSFACTORY_H_

#include "model/ObsStreamQG.h"
#include "model/ObsWindQG.h"
#include "model/ObsWSpeedQG.h"
#include "model/ObservationsQG.h"

namespace qg {

void instantiateQgObsFactory() {
  static ObsMaker<ObsStreamQG> makerStream_("Stream");
  static ObsMaker<ObsWindQG>   makerWind_("Wind");
  static ObsMaker<ObsWSpeedQG> makerWSpeed_("WSpeed");
}

}  // namespace qg

#endif  // QG_MODEL_INSTANTIATEQGOBSFACTORY_H_
