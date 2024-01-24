/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_INSTANTIATEINFLATIONFACTORY_H_
#define OOPS_BASE_INSTANTIATEINFLATIONFACTORY_H_

#include "oops/base/MultiplicativeInflation.h"
#include "oops/base/RTPP.h"
#include "oops/base/RTPS.h"

namespace oops {

template <typename MODEL> void instantiateInflationFactory() {
  static InflationMaker<MODEL, MultiplicativeInflation<MODEL> >  makerMult_("Multiplicative");
  static InflationMaker<MODEL, RTPP<MODEL> >  makerRTPP_("RTPP");
  static InflationMaker<MODEL, RTPS<MODEL> >  makerRTPS_("RTPS");
}

}  // namespace oops

#endif  // OOPS_BASE_INSTANTIATEINFLATIONFACTORY_H_
