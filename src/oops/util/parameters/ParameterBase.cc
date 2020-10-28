/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/util/parameters/ParameterBase.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

ParameterBase::ParameterBase(Parameters *parent) {
  if (parent) {
    parent->registerChild(*this);
  }
}

}  // namespace oops
