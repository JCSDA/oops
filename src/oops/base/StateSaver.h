/*
 * (C) Copyright 2023 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <memory>

#include "oops/base/PostBase.h"
#include "oops/base/Variables.h"

#include "eckit/config/LocalConfiguration.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"

namespace oops {

template <typename FLDS> class StateSaver : public PostBase<FLDS> {
 public:
  explicit StateSaver(const eckit::Configuration & conf)
    : PostBase<FLDS>(conf), vars_(conf.getStringVector("variables")), save_()
  {}
  ~StateSaver() {}

  const FLDS & getState() const {return *save_;}
  FLDS & getState() {return *save_;}

 private:
  void doProcessing(const FLDS & xx) override {save_.reset(new FLDS(vars_, xx));}

  oops::Variables vars_;
  std::unique_ptr<FLDS> save_;
};

}  // namespace oops

