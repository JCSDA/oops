/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <ostream>
#include <string>
#include <vector>

#include "atlas/util/Config.h"

#include "oops/util/Logger.h"

#include "lorenz95/ModelData.h"
#include "lorenz95/Resolution.h"

// -------------------------------------------------------------------------------------------------

namespace lorenz95 {

// -------------------------------------------------------------------------------------------------

ModelData::ModelData(const Resolution & geometry) {}

// -------------------------------------------------------------------------------------------------

ModelData::~ModelData() {}

// -------------------------------------------------------------------------------------------------

const oops::Variables ModelData::defaultVariables() {
    return oops::Variables(std::vector<std::string>({"x"}));
}

// -------------------------------------------------------------------------------------------------

const eckit::LocalConfiguration ModelData::modelData() const {
  return eckit::LocalConfiguration();
}

// -------------------------------------------------------------------------------------------------

void ModelData::print(std::ostream & os) const {
  os << "lorenz95::ModelData::modelData(): " << modelData();
}

// -------------------------------------------------------------------------------------------------

}  // namespace lorenz95
