/*
 * (C) Copyright 2018-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include<memory>
#include<string>
#include<vector>

#include "oops/base/LinearModelBase.h"
#include "oops/generic/HybridLinearModel.h"
#include "oops/generic/IdentityLinearModel.h"
#include "oops/interface/LinearModel.h"

namespace oops {

template <typename MODEL> void instantiateLinearModelFactory() {
  static LinearModelMaker<MODEL, IdentityLinearModel<MODEL> > makerIdentityLinearModel_("Identity");
  static LinearModelMaker<MODEL, HybridLinearModel<MODEL> > makerHybridTangentLinearModel_("HTLM");

  typedef LinearModelMaker<MODEL, interface::LinearModel<MODEL>>   Maker_;
  static std::vector<std::shared_ptr<Maker_>> makers_;
  for (const std::string & name : MODEL::LinearModel::names()) {
    std::shared_ptr<Maker_> pp(new Maker_(name));
    makers_.push_back(pp);
  }
}

}  // namespace oops

