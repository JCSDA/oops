/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_INSTANTIATEQGLOCALIZATIONFACTORY_H_
#define QG_MODEL_INSTANTIATEQGLOCALIZATIONFACTORY_H_

#include "oops/interface/LocalizationBase.h"

#include "oops/qg/LocalizationMatrixQG.h"
#include "oops/qg/QgTraits.h"

namespace qg {

void instantiateQgLocalizationFactory() {
  static oops::interface::LocalizationMaker<qg::QgTraits, LocalizationMatrixQG> makerQG_("QG");
}

}  // namespace qg

#endif  // QG_MODEL_INSTANTIATEQGLOCALIZATIONFACTORY_H_
