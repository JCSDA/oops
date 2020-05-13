/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_GENERIC_INSTANTIATELOCALIZATIONFACTORY_H_
#define OOPS_GENERIC_INSTANTIATELOCALIZATIONFACTORY_H_

#include "oops/base/LocalizationBase.h"
#include "oops/generic/LocalizationID.h"

namespace oops {

template <typename MODEL> void instantiateLocalizationFactory() {
  static oops::LocalizationMaker<MODEL, LocalizationID<MODEL> >  makerID_("ID");
}

}  // namespace oops

#endif  // OOPS_GENERIC_INSTANTIATELOCALIZATIONFACTORY_H_
