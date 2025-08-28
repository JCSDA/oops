/*
 * (C) Copyright 2018-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include<memory>
#include<string>
#include<vector>

#include "oops/base/ModelBase.h"
#include "oops/generic/IdentityModel.h"
#include "oops/generic/PseudoModel.h"
#include "oops/interface/Model.h"

namespace oops {

template <typename MODEL> void instantiateModelFactory() {
  static ModelMaker<MODEL, IdentityModel<MODEL> > makerIdentityModel_("Identity");
  static ModelMaker<MODEL, PseudoModel<MODEL> > makerPseudoModel_("PseudoModel");

  typedef ModelMaker<MODEL, interface::Model<MODEL>>   Maker_;
  static std::vector<std::shared_ptr<Maker_>> makers_;
  for (const std::string & name : MODEL::Model::names()) {
    std::shared_ptr<Maker_> pp(new Maker_(name));
    makers_.push_back(pp);
  }
}

}  // namespace oops

